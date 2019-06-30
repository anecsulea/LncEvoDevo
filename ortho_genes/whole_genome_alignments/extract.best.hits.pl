use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################
##############################################################

sub readAlignmentStats{
    my $pathin=$_[0];
    my $sp1=$_[1];
    my $sp2=$_[2];
    my $hits12=$_[3];
    my $hits21=$_[4];

    open(my $input, $pathin);
    my $line=<$input>; ## header
    chomp $line;
    my @s=split("\t", $line);
    
    my %header;
    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }

    $line=<$input>;
    
    while($line){
	chomp $line;
	my @s=split("\t", $line);
	
	my $id1=$s[$header{"ID.".$sp1}];
	my $id2=$s[$header{"ID.".$sp2}];

	my $lenungapped=$s[$header{"LengthUngapped"}];
	my $lenidentical=$s[$header{"LengthIdentical"}];

	if(exists $hits12->{$id1}){
	    $hits12->{$id1}{$id2}={"ungapped"=>$lenungapped, "identical"=>$lenidentical};
	} else{
	    $hits12->{$id1}={$id2=>{"ungapped"=>$lenungapped, "identical"=>$lenidentical}};
	}

	if(exists $hits21->{$id2}){
	    $hits21->{$id2}{$id1}={"ungapped"=>$lenungapped, "identical"=>$lenidentical};
	} else{
	    $hits21->{$id2}={$id1=>{"ungapped"=>$lenungapped, "identical"=>$lenidentical}};
	}
	
	$line=<$input>;
    }

    close($input);
}

##############################################################

sub readExonBlocks{
    my $pathin=$_[0];
    my $blocks=$_[1];

    open(my $input, $pathin);
    
    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my $gene=$s[0];
	my $chr=$s[2];
	my $start=$s[3];
	my $end=$s[4];
	my $strand=$s[5];

	if(!(exists $blocks->{$gene})){
	    $blocks->{$gene}={"chr"=>$chr, "strand"=>$strand, "start"=>[], "end"=>[]};
	}

	push(@{$blocks->{$gene}{"start"}}, $start);
	push(@{$blocks->{$gene}{"end"}}, $end);
	
	$line=<$input>;
    }

    close($input);
}

##############################################################

sub computeExonicLength{
    my $blocks=$_[0];
    my $exlen=$_[1];

    foreach my $gene (keys %{$blocks}){
	my $totlen=0;
	my $nb=@{$blocks->{$gene}{"start"}};

	for(my $i=0; $i<$nb; $i++){
	    $totlen+=(${$blocks->{$gene}{"end"}}[$i]-${$blocks->{$gene}{"start"}}[$i]+1);
	}

	$exlen->{$gene}=$totlen;
    }
}

##############################################################

sub extractBestHits{
    my $hits=$_[0];
    my $utr1=$_[1];
    my $utr2=$_[2];
    my $exoniclength=$_[3];
    my $minalnfraction=$_[4];
    my $minratiosecondbest=$_[5];
    my $besthits=$_[6];

    foreach my $id1 (keys %{$hits}){
	if(!(exists $utr1->{$id1})){
	    if(!exists $exoniclength->{$id1}){
		print "Weird! cannot find exonic length for ".$id1."\n";
		exit(1);
	    }
	    
	    my $totlen=$exoniclength->{$id1};
	    
	    if($totlen<=0){
		print "Weird! length should be positive for ".$id1." ".$totlen."\n";
		exit(1);
	    }
	    
	    my %alnlengths;
	    
	    foreach my $id2 (keys %{$hits->{$id1}}){
		if(!(exists $utr2->{$id2})){
		    my $thislen=$hits->{$id1}{$id2}{"identical"}; ## length identical sequence
		    
		    if(exists $alnlengths{$thislen}){
			push(@{$alnlengths{$thislen}}, $id2); 
		    } else{
			$alnlengths{$thislen}=[$id2]; 
		    }
		}
	    }

	    my $nbaln=keys %alnlengths;

	    if($nbaln>0){
		my @uniquelengths=keys %alnlengths;
		my @sortedlengths=sort {$a<=>$b} @uniquelengths;
		my $nbunique=@uniquelengths;
		
		my $maxlen=$sortedlengths[-1];
		my $nbmax=@{$alnlengths{$maxlen}};
		
		my $thisalnfraction=($maxlen+0.0)/($totlen+0.0);
		
		if($nbmax==1 && $thisalnfraction>=$minalnfraction && $maxlen>0){
		    my $best=${$alnlengths{$maxlen}}[0];
		    
		    if($nbunique==1){
			$besthits->{$id1}=$best;
		    }  else{
			my $secondbestlength=$sortedlengths[$nbunique-2];
			
			if($secondbestlength==0){
			    $besthits->{$id1}=$best;
			}
			else{
			    my $thisratio=($maxlen+0.0)/($secondbestlength+0.0);
			    
			    if($thisratio>=$minratiosecondbest){
				$besthits->{$id1}=$best;
			    }
			}
		    }
		}
	    }
	}
    }
}

##############################################################

sub readUTR{
    my $pathin=$_[0];
    my $utrs=$_[1];

    open(my $input, $pathin);
    my $line;
    chomp $line;
    my @s=split("\t", $line);
    
    my %header;
    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }

    $line=<$input>;
    
    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $gid=$s[$header{"GeneID"}];

	$utrs->{$gid}=1;
	
	$line=<$input>;
    }

    close($input);
}

##############################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script extract best hits from projection alignments.\n";
    print "\n";
    
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

##############################################################
##############################################################

my %parameters;

$parameters{"species1"}="NA";
$parameters{"species2"}="NA";
$parameters{"pathExonBlocks1"}="NA";
$parameters{"pathExonBlocks2"}="NA";
$parameters{"pathUTR1"}="NA";
$parameters{"pathUTR2"}="NA";
$parameters{"pathAlignmentStats"}="NA";
$parameters{"minAlignedFraction"}="NA";
$parameters{"minRatioSecondBest"}="NA";
$parameters{"pathBestHits12"}="NA";
$parameters{"pathBestHits21"}="NA";
$parameters{"pathReciprocalBestHits"}="NA";

my @defaultpars=("species1", "species2", "pathExonBlocks1", "pathExonBlocks2", "pathUTR1", "pathUTR2", "pathAlignmentStats", "minAlignedFraction", "minRatioSecondBest", "pathBestHits12", "pathBestHits21",  "pathReciprocalBestHits");

my %defaultvalues;

foreach my $par (keys %parameters){
    $defaultvalues{$par}=$parameters{$par};
}

## update arguments

my $nbargs=@ARGV;

for(my $i=0;$i<$nbargs; $i++){
    my $arg=$ARGV[$i];
    $arg=substr $arg,2;
    
    my @s=split("=",$arg);
    my $parname=$s[0];
    my $parval=$s[1];
    
    if(exists $parameters{$parname}){
	$parameters{$parname}=$parval;
    }
    else{
	print "Error: parameter ".$parname." was not recognized!!!\n";
	printHelp(\@defaultpars, \%defaultvalues);
	exit(1);
    }
}

## show parameters

print "\n";

print "Running program with the following parameters:\n";

foreach my $par (@defaultpars){
    print "--".$par."=".$parameters{$par}."\n";
}

print "\n";

##############################################################
##############################################################

my $sp1=$parameters{"species1"};
my $sp2=$parameters{"species2"};

print "First species: ".$sp1.".\n";
print "Second species: ".$sp2.".\n";

##############################################################

print "Reading exon blocks and computing total exonic length...\n";

my %exonblocks1;
readExonBlocks($parameters{"pathExonBlocks1"}, \%exonblocks1);

my %exoniclength1;
computeExonicLength(\%exonblocks1, \%exoniclength1);

my $nbg=keys %exonblocks1;
my $nbel=keys %exoniclength1;

print "Found exon blocks for ".$nbg." genes, first species.\n";
print "Computed exonic length for ".$nbel." genes, first species.\n";


my %exonblocks2;
readExonBlocks($parameters{"pathExonBlocks2"}, \%exonblocks2);

my %exoniclength2;
computeExonicLength(\%exonblocks2, \%exoniclength2);


my $nbg=keys %exonblocks2;
my $nbel=keys %exoniclength2;

print "Found exon blocks for ".$nbg." genes, second species.\n";
print "Computed exonic length for ".$nbel." genes, second species.\n";

print "Done.\n";

##############################################################

print "Reading alignment statistics...\n";

my %hits12;
my %hits21;

readAlignmentStats($parameters{"pathAlignmentStats"}, $sp1, $sp2, \%hits12, \%hits21);
   
print "Done.\n";

##############################################################

print "Reading putative UTRs...\n";
my %utr1;
readUTR($parameters{"pathUTR1"}, \%utr1);
my $nbu1=keys %utr1;

print "Found ".$nbu1." putative UTRs for ".$sp1."\n";

my %utr2;
readUTR($parameters{"pathUTR2"}, \%utr2);

my $nbu2=keys %utr2;

print "Found ".$nbu2." putative UTRs for ".$sp2."\n";

print "Done.\n";

##############################################################

print "Extracting best hits...\n";

my $minratio=$parameters{"minRatioSecondBest"}+0.0;

if($minratio<=1){
    print "Weird! ratio best-to-second-best <=1!\n";
    exit(1);
}

print "Minimum ratio of ungapped length best to second best: ".$minratio."\n";

my %besthits12;
my %besthits21;

my $minalnfraction=$parameters{"minAlignedFraction"}+0.0;
print "Minimum aligned fraction: ".$minalnfraction."\n";

extractBestHits(\%hits12, \%utr1, \%utr2, \%exoniclength1, $minalnfraction, $minratio, \%besthits12);
extractBestHits(\%hits21, \%utr2, \%utr1, \%exoniclength2, $minalnfraction, $minratio, \%besthits21);

my $nbh12=keys %besthits12;
my $nbh21=keys %besthits21;

print "Found ".$nbh12." best hits for ".$sp1." to ".$sp2."\n";
print "Found ".$nbh21." best hits for ".$sp2." to ".$sp1."\n";

print "Done.\n";

##############################################################

print "Writing output for best hits...\n";

open(my $output, ">".$parameters{"pathBestHits12"});
print $output "ID.".$sp1."\tID.".$sp2."\n";

foreach my $id1 (keys %besthits12){
    print $output $id1."\t".$besthits12{$id1}."\n";
}

close($output);

open(my $output, ">".$parameters{"pathBestHits21"});
print $output "ID.".$sp2."\tID.".$sp1."\n";

foreach my $id2 (keys %besthits21){
    print $output $id2."\t".$besthits21{$id2}."\n";
}

close($output);
print "Done.\n";

##############################################################

print "Extracting reciprocal best hits and writing output...\n";

open(my $output, ">".$parameters{"pathReciprocalBestHits"});

print $output "ID.".$sp1."\tID.".$sp2."\n";

foreach my $id1 (keys %besthits12){
    my $id2=$besthits12{$id1};

    if(exists $besthits21{$id2} && $besthits21{$id2} eq $id1){
	print $output $id1."\t".$id2."\n";
    }
}

close($output);

print "Done.\n";

##############################################################


