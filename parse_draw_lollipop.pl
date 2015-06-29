use strict;

#my @bams = </mnt/home2/zhiyi/APOBEC/E14_APOBEC_SEQ_May2015/051515_pilot_run/control_analysis/bam/lambda/*.bam>;
my @bams1 = </mnt/home2/zhiyi/APOBEC/E14_APOBEC_SEQ_May2015/051515_pilot_run/control_analysis/bam/*/*.bam>;
my @bams2 =</mnt/home2/zhiyi/APOBEC/E14_APOBEC_SEQ_May2015/052915_tech_replicate_run/control_analysis/bam/*bam>;
my @bams = (@bams1,@bams2);

foreach my $bam (@bams)
{
    my $generic = $bam;
    $generic =~ s/.*\///g;
    $generic =~ s/\.bam//;
    print STDERR "dealing with $generic\n";
    my $sam = $generic.".sam";

    my $c1 = "samtools view -f 98 $bam > $sam";
    my $out = "unconverted_".$generic.".sam";
    system($c1);
    my $c2 = "perl convert_to_biq_analyzer.pl $sam";
    system($c2);
    
}
    
