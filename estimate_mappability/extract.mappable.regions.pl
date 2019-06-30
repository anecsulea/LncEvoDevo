use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################

sub readCoordinates{
     my $path=$_[0];
     my $coords=$_[1];

     open(my $input, "zcat $path |");
     my $line=<$input>;
     my $nbread=0;
     
     while($line){
	 my $prefix=substr $line,0,1;
	 
	 if($prefix eq ">"){
	     $line=substr $line,1;
	     chomp $line;
	     my @s=split(" ",$line);
	     my $id=$s[0];
	     my $c=$s[1];
	     my @t=split(":",$c);
	   
	     my $chr=$t[0];

	     my @p=split("-",$t[1]);
	     my $start=$p[0]+1; ## coordinates now start at 1, like in the sam file

	     $coords->{$id}={"chr"=>$chr,"start"=>$start};

	     $nbread++;
	     if($nbread%1000000==0){
		 print "read ".$nbread." fake reads.\n";
	     }
	 }

	 $line=<$input>;
     }

     close($input);
}

##############################################################

sub checkAlignments{
    my $pathin=$_[0];
    my $coords=$_[1];
    my $readlen=$_[2];
    
    my @s=split("\\.",$pathin);
    my $ext=$s[-1];
    my $input;
    
    if($ext eq "bam"){
	open($input, "samtools view $pathin |");
    }
    else{
	if($ext eq "sam"){
	    open($input, $pathin);
	}
	else{
	    if($ext eq "gz"){
		open($input, "zcat $pathin |");	
	    }
	    else{
		print "unknown file extension for ".$pathin."\n";
		exit(1);
	    }
	}
    }

    my $line=<$input>;

    my $nbwrongchr=0;
    my $nbwrongpos=0;
    my $nbwrongstrand=0;
    my $nbambiguous=0;
    my $nblocalaln=0;
    my $nbok=0;
    my $nbunmapped=0;
    my $nbdone=0;
    
    my %removedreads;

    my $okcigar=$readlen."M"; ### no indels, no clipping - stringent settings

    my $firstchar=substr $line,0,1;

    while($firstchar eq "@"){
	$line=<$input>;
	$firstchar=substr $line,0,1;
    }

    while($line){
	my $thisreadok=0;

	$nbdone++;
	if($nbdone%1000000==0){
	    print $nbdone." alignments read.\n";
	}
	chomp $line;
	my @s=split("\t",$line);
	my $readid=$s[0];
	my $flag=$s[1]+0;
	my $chr=$s[2];
	my $start=$s[3]+0;
	my $cigar=$s[5];
	my @AS=grep(/AS:i:/,@s);
	my @XS=grep(/XS:i:/,@s);

	my $thisscore=-100;
	my $nextscore=-100;

	if(@AS==1){
	    my @u=split(":",$AS[0]);
	    $thisscore=$u[2]+0;
	}
	else{
	    if($chr ne "*"){
		print "saw weird AS:i ".join("\t",@AS)." in ".$line."\n";
		exit(1);
	    }
	}

	if(@XS==1){
	    my @u=split(":",$XS[0]);
	    $nextscore=$u[2]+0;
	}

	
	if(exists $coords->{$readid}){
	    my $realstart=$coords->{$readid}{"start"};
	    my $realchr=$coords->{$readid}{"chr"};
	    
	    if($chr eq $realchr){
		if($realstart==$start){
		    if(!($flag & 16)){
		
			if($thisscore>$nextscore){
			    if($cigar eq $okcigar){
				$thisreadok=1;
				$nbok++;
			    }
			    else{
				$nblocalaln++;
			    }
			}
			else{
			    $nbambiguous++;
			}
		    }
		    else{
			$nbwrongstrand++;
		    }
		}
		else{
		    $nbwrongpos++;
		}
	    }
	    else{
		if($chr eq "*"){
		    $nbunmapped++;  
		}
		else{
		    $nbwrongchr++;
		}
	    }
	    
	    if($thisreadok==0){
		delete $coords->{$readid};
		$removedreads{$readid}=1;
	    }
	}
	else{
	    if(!(exists $removedreads{$readid})){
		print "Weird! cannot find ".$readid." in the original fake reads data, and we didn't remove it previously.\n";
		exit(1);
	    }
	}

	$line=<$input>;
    }

    close($input);

    print "There were ".$nbok." correct alignments.\n";
    print "There were ".$nbunmapped." unmapped reads (likely containing Ns).\n";
    print "There were ".$nbwrongchr." reads mapped on the wrong chromosome.\n";
    print "There were ".$nbwrongpos." reads mapped at the wrong position.\n";
    print "There were ".$nbwrongstrand." reads mapped on the wrong strand.\n";
    print "There were ".$nbambiguous." reads mapped ambiguously.\n";
    print "There were ".$nblocalaln." reads without end-to-end alignments.\n";
    
}

##############################################################

sub extractMappedRegions{
    my $coords=$_[0];
    my $nbtotreads=$_[1];
    my $readlen=$_[2];
    my $regions=$_[3];
 
    my $currentchr="NA";
    my $currentstart="NA";
    my $currentend="NA";

    print "Total number of reads ".$nbtotreads."\n";
   
    for(my $i=0; $i<$nbtotreads; $i++){	## this is very important because the reads are ordered !!!
	my $id="read".$i;

	if(exists $coords->{$id}){
	    my $thischr=$coords->{$id}{"chr"};
	    my $thisstart=$coords->{$id}{"start"};
	    my $thisend=$thisstart+$readlen-1;

	    if($currentend eq "NA" || $currentstart eq "NA"){
		$currentstart=$thisstart;
		$currentend=$thisend;
		$currentchr=$thischr;
	    }
	    else{
		if($currentchr eq $thischr && $thisstart>=$currentstart && $thisstart<=($currentend+1)){
		    $currentend=$thisend;
		}
		else{
		    if(exists $regions->{$currentchr}){
			push(@{$regions->{$currentchr}{"start"}},$currentstart);
			push(@{$regions->{$currentchr}{"end"}},$currentend);
		    }
		    else{
			$regions->{$currentchr}={"start"=>[$currentstart],"end"=>[$currentend]};
		    }

		    $currentstart=$thisstart;
		    $currentend=$thisend;
		    $currentchr=$thischr;
		}
	    }
	}
    }
    
    if(exists $regions->{$currentchr}){
	push(@{$regions->{$currentchr}{"start"}},$currentstart);
	push(@{$regions->{$currentchr}{"end"}},$currentend);
    }
    else{
	$regions->{$currentchr}={"start"=>[$currentstart],"end"=>[$currentend]};
    }
}

##############################################################

sub printHelp{
    
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script extracts mappable regions of the genome.\n";
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
$parameters{"pathReads"}="NA";
$parameters{"pathAlignment"}="NA";
$parameters{"readLength"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathReads","pathAlignment","readLength","pathOutput");
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

print "Reading original coordinates...\n";
my %coordinates;
readCoordinates($parameters{"pathReads"},\%coordinates);
my $nbtotreads=keys %coordinates;
print "There are ".$nbtotreads." reads in total.\n";
print "Done.\n";

print "Reading alignment and checking coordinates...\n";

my $readlen=$parameters{"readLength"}+0;
print "Read length ".$readlen."\n";
checkAlignments($parameters{"pathAlignment"},\%coordinates, $readlen);
my $nbkept=keys %coordinates;
print "There are ".$nbkept." reads with correct alignments.\n";
print "Done.\n";

##############################################################

print "Extracting mapped regions...\n";
my %mapped;
extractMappedRegions(\%coordinates, $nbtotreads, $readlen, \%mapped);
print "Done.\n";

print "Writing output...\n";
open(my $output, ">".$parameters{"pathOutput"});
foreach my $chr (keys %mapped){
    my $nbreg=@{$mapped{$chr}{"start"}};
    
    for(my $i=0; $i<$nbreg; $i++){
	my $start=${$mapped{$chr}{"start"}}[$i];
	my $end=${$mapped{$chr}{"end"}}[$i];
	print $output $chr."\t".$start."\t".$end."\n";
    }
}
close($output);
print "Done.\n";

##############################################################
