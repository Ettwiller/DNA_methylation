#!/usr/bin/perl


use strict;
use warnings;
use Getopt::Long qw(GetOptions);


my $samfile;
my $path_to_BiQ = "/home/ettwiller/exe/BiQ_analyser_diagram/MethylationDiagrams.jar";
my $random = 0;

my $usage_sentence = "perl $0 --sam samfile.sam OPTIONAL : --random 1\n";

GetOptions ("sam=s" => \$samfile,    # numeric
	    "random=s" => $random
    ) or die $usage_sentence;

if (!$samfile) {die $usage_sentence;}

my $generic = $samfile;
$generic =~ s/.*\///g;
$generic =~ s/\.sam//;
my $sam = "tmp_".$generic.".sam";

my $c1 = "samtools view -f 98 $samfile > $sam";
system($c1);


#create a random sam for comparison
my ($methylation, $random_methylation) = random_sam($sam);
unlink($sam);
my ($biq_methylation) = biq_format($methylation);
my ($biq_methylation_random) = biq_format($random_methylation);

my $generic = $sam;
$generic =~ s/.*\///g; 
$generic =~ s/\.sam//;
#$generic = $generic."_random"; 
my $IN = $generic.".txt";

open (IN, ">$IN") or die;
my $size = @$biq_methylation;
print IN ">testing\n";
my $cutoff = 0;
for (my $i=0; $i<$size; $i++)
{
    my $line = $$biq_methylation[$i];
#    my $line = $$biq_methylation_random[$i];
 
    my $pattern = "I"; #get the converted 
    my $total_pattern = "[I, 0]";
    my $c = () = $line =~ /$pattern/g;
    my $C= () = $line =~ /$total_pattern/g;
    if ($c>0 && $cutoff < 200)
    {
#	$line = $random;
	#inverse the colors !!!!
#	$line=~ s/I/z/g; 
#	$line=~ s/0/I/g;
#        $line=~ s/z/0/g;
	$cutoff++;
	print IN "Seq $i\t$line\n";
    }
   
}

my $txt = $generic.".txt".".tmp";
my $html = $generic.".html".".tmp";
my $png1 = $generic.".png";
my $svg1 = $generic.".svg";
my $png2 = $generic.".png".".tmp";
my $svg2 = $generic.".svg".".tmp";

my $command  = "java -jar /home/ettwiller/exe/BiQ_analyser_diagram/MethylationDiagrams.jar $IN uncompressed $txt $html $png1 $svg1 $png2 $svg2 true false 1 5000 5 false true FFFF00 0000FF 808080";

system($command);
my $command2 = "rm $generic*tmp";
system($command2);

sub biq_format{
    my ($array) =@_;
    my @tmp; my @biq; my $MAX=0;
    foreach my $line(@$array){
	
	$line=~ s/[x,h,z]/I/g; #I = converted
	$line=~ s/[X,H,Z]/0/g; #0 = unconverted
	$line=~ s/\./\-/g;
	my $size = length($line);
	if ($size > $MAX){ $MAX=$size; }
	push @tmp, $line;
	
    }
    foreach my $line (@tmp)
    {
	my $add = $MAX - length($line);
	my $addition = "-" x $add;
	$line  = $line.$addition;
	push @biq, $line;
    }

    return (\@biq);
}
    



sub random_sam {
    my ($sam)=@_;
    my @real; my @random;

    open (SAM, $sam) or die "can't open $sam file\n";
    my $raw;
    foreach my $line (<SAM>)
    {
	chomp $line;
	my @tmp = split /\t/, $line;
	my $methylation = $tmp[13];
	$methylation =~ s/XM:Z://;
	my $result = get_raw_number($methylation);
	foreach my $pattern(keys %$result){$$raw{$pattern} += $$result{$pattern};}
    }
    close SAM;
    open (SAM, $sam) or die "can't open $sam file\n";
    foreach my $line (<SAM>)
    {
	chomp $line;
	my @tmp = split /\t/, $line;
	my $methylation = $tmp[13];
	$methylation =~ s/XM:Z://;
	my $random_methylation = generage_random_methylation($methylation, $raw);
	push @random, $random_methylation;
	push @real, $methylation;
    }
    return(\@real, \@random);
}


sub generage_random_methylation {
    my ($line, $raw)=@_;
    my $randomized_methylation;
    my @tmp=split//, $line;
    foreach my $e (@tmp)
    {

	if ($e =~ /[a-z,A-Z]/)
	{
	    my $UC = uc($e);
	    my $LC = lc($e);
	  
	    my $ratio = $$raw{$UC}/($$raw{$LC} + $$raw{$UC});
	    my $rand = rand(1);
	    if ($rand > $ratio){ $randomized_methylation = $randomized_methylation.$LC;}
	    else {$randomized_methylation = $randomized_methylation.$UC;}
				 
	}
	else {
	    $randomized_methylation = $randomized_methylation.$e;
	}
    }
    return $randomized_methylation;
}


sub get_raw_number {
    my ($line) =@_;
    my %result;
    my @all = ("h", "H", "x", "X", "z", "Z");
    foreach my $pattern (@all)
    {
	my $c = () = $line =~ /$pattern/g;
	$result{$pattern}=$c;
    }
    return (\%result);
}
