use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##################################################################

sub readGTF{
    my $pathin=$_[0];
    my $exoncoords=$_[1];
    my $genetx=$_[2];
    my $txex=$_[3];
 
    open(my $input, $pathin);
    
    my $line=<$input>;
    my $prefix=substr $line, 0,1;
    
    while($prefix eq "#"){
	$line=<$input>;
	$prefix=substr $line, 0,1;
    }

    my $nbdiscarded=0;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	
	my $type=$s[2];

	if($type eq "exon"){
	
	    my $chr=$s[0];
	    my $start=$s[3]+0; ## 1-based
	    my $end=$s[4]+0;
	    my $strand=$s[6];

	    if($strand eq "+"){
		$strand="1";
	    } else{
		if($strand eq "-"){
		    $strand="-1";
		} else{
		    if($strand ne "."){
			print "Unknown strand ".$strand." at line ".$line."\n";
			exit(1);
		    }
		}
	    }

	    if($strand ne "."){
		my $info=$s[8];
		my @t=split(";", $info);
		my $geneid=findInfo("gene_id", \@t);
		my $txid=findInfo("transcript_id", \@t);
		my $classcode=findInfo("class_code", \@t);
		my $oldid=findInfo("oId", \@t);
	
		if($txid eq "NA"){
		    print "could not find transcript in ".$line."\n";
		    exit(1);
		}		

		if($geneid eq "NA"){
		    print "could not find gene in ".$line."\n";
		    exit(1);
		}
		
		my $exonid=$chr.",".$start.",".$end.",".$strand;
		
		## fill in exon coords
		
		$exoncoords->{$exonid}={"chr"=>$chr, "start"=>$start, "end"=>$end, "strand"=>$strand};
			
		if(exists $genetx->{$geneid}){
		    $genetx->{$geneid}{$txid}=1;
		}
		else{
		    $genetx->{$geneid}={$txid=>1};
		}
		
		if(exists $txex->{$txid}){
		    $txex->{$txid}{$exonid}=1;
		}
		else{
		    $txex->{$txid}={$exonid=>1};
		}
	    }
	    else{
		$nbdiscarded++; 
	    }
	}
   	
	$line=<$input>;
    }
    
    close($input);
      
    print "We discarded ".$nbdiscarded." exons with undefined strand.\n";
}


##############################################################

sub findInfo{
    my $pattern=$_[0];
    my $array=$_[1];
    
    my $res="NA";
    
    my @grepres=grep(/${pattern}/,@{$array});

    my $nbg=@grepres;
    my $nbreal=0;

    if(@grepres==1){
	my @t=split("\"",$grepres[0]);
	$res=$t[1];
    } else{
	my $nbreal=0;
	
	foreach my $g (@grepres){
	    $g =~ s/^\s+//; ## remove whitespace
	    my @u=split(" ",$g);

	    if($u[0] eq $pattern){
		$nbreal++;
		my @t=split("\"",$g);
		$res=$t[1];
	    }
	}
    }
    
    if($nbreal>1){
	return "NA";
    }
    return $res;
}


##################################################################

sub extractSpliceSites{
    my $genetx=$_[0];
    my $txexons=$_[1];
    my $exoncoords=$_[2];
    my $splice5=$_[3];
    my $splice3=$_[4];

    foreach my $geneid (keys %{$genetx}){
	foreach my $txid (keys %{$genetx->{$geneid}}){
	    my $nbexons=keys %{$txexons->{$txid}};
	    
	    if($nbexons>=2){
		my %hashcoords;
		
		my $chr="NA";
		my $strand="NA";
	    	
		foreach my $exonid (keys %{$txexons->{$txid}}){
		    my $start=$exoncoords->{$exonid}{"start"};
		    my $end=$exoncoords->{$exonid}{"end"};
		    
		    if(exists $hashcoords{$start}){
			print "Weird! already saw exon start ".$start." for ".$txid."\n";
			exit(1);
		    }
		    
		    $hashcoords{$start}=$end;
		    
		    if($chr eq "NA"){
			$chr=$exoncoords->{$exonid}{"chr"};
			$strand=$exoncoords->{$exonid}{"strand"};
		    }
		}
			    
		my @uniquestart=keys %hashcoords;
		my @startex=sort {$a<=>$b} @uniquestart;
		my @endex;
		
		foreach my $start (@startex){
		    push(@endex, $hashcoords{$start});
		}
		
		for(my $i=0; $i<($nbexons-1); $i++){
		    my $start1=$startex[$i];
		    my $end1=$endex[$i];
		    
		    my $start2=$startex[$i+1];
		    my $end2=$endex[$i+1];

		    my $id5="NA";
		    my $id3="NA";
		    
		    if($strand eq "1"){
			$id5=$chr.",".($end1+1).",".($end1+2).",".$strand;
			$id3=$chr.",".($start2-2).",".($start2-1).",".$strand;
		    } else{
			if($strand eq "-1"){
			    $id3=$chr.",".($end1+1).",".($end1+2).",".$strand;
			    $id5=$chr.",".($start2-2).",".($start2-1).",".$strand;
			} else{
			    print "Weird strand: ".$strand." geneid ".$geneid."\n";
			    exit(1);
			}
		    }
		    
		    if(exists $splice5->{$geneid}){
			$splice5->{$geneid}{$id5}=1;
		    } else{
			$splice5->{$geneid}={$id5=>1};
		    }
		    
		    if(exists $splice3->{$geneid}){
			$splice3->{$geneid}{$id3}=1;
		    } else{
			$splice3->{$geneid}={$id3=>1};
		    }
		}
	    }
	}
    }
}

##################################################################

sub readPhastCons{
    my $pathin=$_[0];
    my $refphast=$_[1];
    
    my @s=split("\\.",$pathin);
    my $nbs=@s;
    my $ext=$s[$nbs-1];

    my $input;
    
    if($ext eq "gz"){
	open($input,"zcat $pathin |");
    }
    else{
	open($input,$pathin);
    }
    
    my $line=<$input>;
    my $currentpos="NA";
    
    while($line){
	my $first=substr $line,0,1;
	
	if($first eq "f"){
	    my @s=split(" ",$line);
	    my $start=$s[2];
	    my @t=split("=",$start);
	    $currentpos=$t[1]+0;
	    $line=<$input>;
	    next;
	}
	else{
	    chomp $line;
	    my $val=$line+0.0;

	    $refphast->{$currentpos}=$val;
	    
	    $currentpos++;
	    $line=<$input>;
	}
    }
    
    close($input);
}


##############################################################

sub computePhastConsScore{
    my $splicecoords=$_[0];
    my $refphast=$_[1];
    my $chr=$_[2];
    my $scores=$_[3];
   
    foreach my $gene (keys %{$splicecoords}){
	my $totbases=0;
	my $coveredbases=0;
	my $sumscore=0; 
	my $thischr="NA";
	
	foreach my $id (keys %{$splicecoords->{$gene}}){
	    my @s=split(",", $id);
	    $thischr=$s[0];
	    my $start=$s[1]+0; ## everything is 1-based
	    my $end=$s[2]+0;

	    for(my $i=$start; $i<=$end; $i++){
		$totbases++;
		
		if(exists $refphast->{$i}){
		    $coveredbases++;
		    $sumscore+=$refphast->{$i};
		}
	    }
	}
	
	if($thischr eq $chr){
	    $scores->{$gene}={"totbases"=>$totbases, "coveredbases"=>$coveredbases, "sumscore"=>$sumscore};
	}
    }
}

##################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script computes phastcons scores for splice sites. \n";
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

$parameters{"pathGTF"}="NA";
$parameters{"pathPhastCons"}="NA";
$parameters{"chr"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("pathGTF", "pathPhastCons", "chr", "pathOutput");

my %numericpars;

foreach my $par (keys %parameters){
    $defaultvalues{$par}=$parameters{$par};
}

## check if help was asked 

foreach my $arg (@ARGV){
    if($arg eq "--help"){
	printHelp(\@defaultpars, \%defaultvalues);
	exit(0);
    }
}

## check new parameters

my $nbargs=@ARGV;

for(my $i=0; $i<$nbargs; $i++){
    my $arg=$ARGV[$i];
    $arg=substr $arg,2;
    my @s=split("=",$arg);
    my $parname=$s[0];
    my $parval=$s[1];
    
    if(exists $parameters{$parname}){
	$parameters{$parname}=$parval;
	
	if(exists $numericpars{$parname}){
	    $parameters{$parname}=$parval+0.0;
	}
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

print "Reading annotations...\n";

my %exoncoords;
my %genetx;
my %txexons;

readGTF($parameters{"pathGTF"}, \%exoncoords, \%genetx, \%txexons);

my $nbg=keys %genetx;
my $nbt=keys %txexons;
my $nbe=keys %exoncoords;

print "Found ".$nbg." genes, ".$nbt." transcripts and ".$nbe." exons in the annotations.\n";

print "Done.\n";

##############################################################

print "Extracting splice site coordinates...\n";

my %splice5;
my %splice3;

extractSpliceSites(\%genetx, \%txexons, \%exoncoords, \%splice5, \%splice3);

print "Done.\n";

##############################################################

print "Reading phastcons scores and computing average for splice sites...\n";

my $chr=$parameters{"chr"};
my $path=$parameters{"pathPhastCons"};

print "Chromosome ".$chr."\n";

if(-e $path){
    my %phastCons;
    readPhastCons($path,\%phastCons);
    
    my %scores5;
    my %scores3;
    
    computePhastConsScore(\%splice5, \%phastCons, $chr, \%scores5);
    computePhastConsScore(\%splice3, \%phastCons, $chr, \%scores3);
    
    open(my $output, ">".$parameters{"pathOutput"});

    print $output "Gene\tTotBases5\tTotBases3\tCoveredBases5\tCoveredBases3\tSumScore5\tSumScore3\n";
    
    foreach my $gene (keys %genetx){
	if(exists $scores5{$gene} && exists $scores3{$gene}){
	    print $output $gene."\t".$scores5{$gene}{"totbases"}."\t".$scores3{$gene}{"totbases"}."\t".$scores5{$gene}{"coveredbases"}."\t".$scores3{$gene}{"coveredbases"}."\t".$scores5{$gene}{"sumscore"}."\t".$scores3{$gene}{"sumscore"}."\n";
	}
	else{
	    if(exists $scores5{$gene} || exists $scores3{$gene}){
		print "Weird ! found one but not both scores for ".$gene."\n";
		exit(1);
	    }
	}
    }
    
    close($output);
}

print "Done.\n";

##############################################################
