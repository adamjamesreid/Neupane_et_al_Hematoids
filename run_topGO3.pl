#!/usr/bin/perl -w
# run_topGO3.pl now outputs the up/down gene ids associated with each GO term
# and runs all three types of GO at once

use strict;
use Getopt::Long;

my $p_cut = 0.05;
my $node_size = 5;
my $method = 'weight01'; # weight01, elim or classic

my $r_prog = '/usr/bin/R';

my $help = 0;
my $clust;
my $go;
my $output_prefix;
my $ann;
my $keep = 0;

#my @types = ('bp', 'mf', 'cc');
my @types = ('bp', 'mf', 'cc');

GetOptions
(
        "h|help"        => \$help,
        "l|list:s"      => \$clust,
        "g|go_terms:s"  => \$go,
        "o|output:s"    => \$output_prefix,
        "p:f"           => \$p_cut,
        "n|node_size:i" => \$node_size,
        "m|method:s"    => \$method,
        "k|keep"        => \$keep,
        "a|ann:s"               => \$ann
);

if(!defined $output_prefix)
{
        $output_prefix = 'topgo_out';
}

if(!defined $go || !defined $clust)
{
        print_usage();
        exit;
}

my %ann = ();

if(defined $ann)
{
        open(A, "<$ann") or die "$!";

        while(<A>)
        {
                chomp;

                my($id, $val) = split /\t/;

                $ann{$id} = $val;
        }
        close A;
}

open(CL, "<$clust") or die "$!";

my @vecs = ();
my @ids = ();
my %cluster = ();

while(<CL>)
{
        chomp;

        push @ids, $_;
        $cluster{$_} = 1;
}
close CL;

my $ids = join "\",\"", @ids;
$ids = '"'.$ids.'"';
push @vecs, $ids;

my $r_file = "$output_prefix\.R";

r_code(\@vecs, $output_prefix, $r_file, $method);

`$r_prog --no-save < $r_file`;

unless($keep)
{
        system("rm -f $r_file") == 0 or die "$!";
}
my @header = ();
my %data = ();

open(OUT, ">$output_prefix\_final.dat") or die "$!";

foreach my $type (@types)
{
        open(FILE, "<$output_prefix\_$type\.dat") or die "$!: $output_prefix\_$type\.dat";

        my %ids = ();
        my @term_order = ();

        while(<FILE>)
        {
                chomp;

                if(/\[\d+\] (.*)/)
                {
                        #my(@ids) = split /\t/, $1;

                        while(/\"([\S]+)\"/g)
                        {
                                $ids{$1} = 1;
                        }

                }
                elsif(/^"\d/)
                {
                        my($row, $id, $term, $ann, $sig, $exp, $p) = split /\t/;

                        $p =~ s/\"//g;
                        $id =~ s/\"//g;
                        $term =~ s/\"//g;

                        next unless $p <= $p_cut;

                        #print "$id\t$term\t$ann\t$sig\t$exp\t$p\t$gene_list\n";
                        $data{$id} = "$id\t$term\t$ann\t$sig\t$exp\t$p";
                        push @term_order, $id;
                }
        }
        close FILE;

        unless($keep)
        {
                system("rm -f $output_prefix\_$type\.dat") == 0 or die "$!";
        }

	open(GENES, "<$output_prefix\_terms_$type\.dat") or die "$!";

        my $term = '';
        my $flag = 0;
        my %list = ();

        while(<GENES>)
        {
                chomp;

                if(/^\$\`(GO:\d+)/)
                {
                        $term = $1;

                        if(exists $data{$term})
                        {
                                $flag = 1;
                        }
                        else
                        {
                                $flag = 0;
                        }
                }
                elsif(/\[/ && $flag == 1)
                {
                        while(/\"([\S]+)\"/g)
                        {
                                push @{$list{$term}}, $1;

                        }
                }
        }
        close GENES;

        unless($keep)
        {
                system("rm -f $output_prefix\_terms_$type\.dat") == 0 or die "$!: $output_prefix\_terms_$type\.dat";
        }

	foreach my $term (@term_order)
        {
                my @gene_list = ();

                foreach (@{$list{$term}})
                {
                        if(exists $cluster{$_})
                        {
                                if(exists $ann{$_})
                                {
                                        push @gene_list, $ann{$_};
                                }
                                else
                                {
                                        push @gene_list, $_;
                                }
                        }
                }

                my $list = join " ", @gene_list;

                print OUT "$type\t$data{$term}\t$list\n";
        }

}
close OUT;

sub r_code
{
        my ($vec, $output_prefix, $r_code, $method) = @_;

        my @files = ();

        open(OUT, ">$r_code") or die "$!: $r_code";

        print OUT <<RCODE;

        #library(topGO, lib.loc="~ar11/R/library")
        library(topGO)

        ref=read.table (file="$go", stringsAsFactors=FALSE)
        names(ref) = c('id', 'go')
        ref.vec = strsplit(ref\$go, split=',', fixed=T)
        names(ref.vec) <- ref\$id
        all.ids <- ref\$id

RCODE

        foreach my $row (@$vec)
        {
                print OUT <<RCODE;
                vec<-c($row)

                scores <- rep(0, nrow(ref)) # list of scores
                names(scores) <- ref\$id
                scores[ref\$id \%in\% vec] <- 1

                geneSelectionFun <- function(score){
                  return(score >= 1)
                }

                GOdataBP <- new("topGOdata",  ontology = 'BP', allGenes = scores, annot = annFUN.gene2GO, gene2GO = ref.vec, geneSelectionFun = geneSelectionFun, nodeSize = $node_size, description = '')
                GOdataMF <- new("topGOdata",  ontology = 'MF', allGenes = scores, annot = annFUN.gene2GO, gene2GO = ref.vec, geneSelectionFun = geneSelectionFun, nodeSize = $node_size, description = '')
                GOdataCC <- new("topGOdata",  ontology = 'CC', allGenes = scores, annot = annFUN.gene2GO, gene2GO = ref.vec, geneSelectionFun = geneSelectionFun, nodeSize = $node_size, description = '')

                resultTopgoBP <- runTest(GOdataBP,algorithm="$method",statistic="Fisher")
                resultTopgoMF <- runTest(GOdataMF,algorithm="$method",statistic="Fisher")
                resultTopgoCC <- runTest(GOdataCC,algorithm="$method",statistic="Fisher")

                resBP<-GenTable( GOdataBP, topGO = resultTopgoBP, orderBy = "topGO", ranksOf = "fisher", topNodes = 50)
                resMF<-GenTable( GOdataMF, topGO = resultTopgoMF, orderBy = "topGO", ranksOf = "fisher", topNodes = 50)
                resCC<-GenTable( GOdataCC, topGO = resultTopgoCC, orderBy = "topGO", ranksOf = "fisher", topNodes = 50)

                sink("$output_prefix\_terms_bp.dat")
                genesInTerm(GOdataBP)
                sink("$output_prefix\_terms_mf.dat")
                genesInTerm(GOdataMF)
                sink("$output_prefix\_terms_cc.dat")
                genesInTerm(GOdataCC)

                sink("$output_prefix\_bp.dat")
                write.table(resBP, sep = "\t")
                sink("$output_prefix\_mf.dat")
                write.table(resMF, sep = "\t")
                sink("$output_prefix\_cc.dat")
                write.table(resCC, sep = "\t")

RCODE
        }
        close OUT;
}
sub print_usage
{
        print STDERR <<USAGE;

        run_topGO3.pl [options]

        Run topGO and produce BP/MF/CC GO terms enriched in list of genes

        -h | -help              Print this message
        -l | -list              Gene list (required)
        -g | -go_terms          All GO terms for genome (required)
        -a | -ann               File of gene annotation e.g. gene names for easier reading
        -o | -output            Output prefix [topgo_out]
        -p                      Adjusted p-value cutoff [0.05]
        -n | -node_size         Node size [5]
        -m | -method            Method (classic, elim, weight01) [weight01]
        -k | -keep              Keep all intermediate files inc. R script [FALSE]

USAGE
}

