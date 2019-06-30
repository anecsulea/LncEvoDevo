use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################

sub readGTF{
    my $pathin=$_[0];
    my $exons=$_[1];
    
    open(my $input, $pathin);
    
    my $line=<$input>;
    
    while($line){
	chomp $line;
	my @s=split("\t", $line);
	
	my $type=$s[2];

	if($type eq "exon"){
	
	    my $chr=$s[0];
	    my $start=$s[3]+0; ## 1-based
	    my $end=$s[4]+0;
	    my $strand=$s[6];
	   
	    my $exonid=$chr.",".$start.",".$end.",".$strand;

	    my $info=$s[8];
	    my @t=split(";", $info);
	    my $geneid=findInfo("gene_id", \@t);

	    if(exists $exons->{$geneid}){
		$exons->{$geneid}{$exonid}={"chr"=>$chr, "start"=>$start, "end"=>$end, "strand"=>$strand};
	    }
	    else{
		$exons->{$geneid}={$exonid=>{"chr"=>$chr, "start"=>$start, "end"=>$end, "strand"=>$strand}};
	    }
	}

	$line=<$input>;
    }

    close($input);
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

##############################################################

sub makeExonBlocks{

    my $exons=$_[0];
    my $collapse=$_[1];
    my $exonblocks=$_[2];

    foreach my $gene (keys %{$exons}){
	
	my %hashcoords;
	my $chr="NA";
	my $strand="NA";

	foreach my $ex (keys %{$exons->{$gene}}){
	    if($chr eq "NA"){
		$chr=$exons->{$gene}{$ex}{"chr"};
	    }

	    if($strand eq "NA"){
		$strand=$exons->{$gene}{$ex}{"strand"};
	    }
	    
	    my $b=$exons->{$gene}{$ex}{"start"};
	    my $e=$exons->{$gene}{$ex}{"end"};
	    
	    if(exists $hashcoords{$b}){
		if($e>$hashcoords{$b}){
		    $hashcoords{$b}=$e;   
		}
	    }
	    else{
		$hashcoords{$b}=$e;
	    }
	}	    
	    
	if($chr eq "NA" || $strand eq "NA"){
	    print $gene." doesn't have any exons, this is weird!!\n";
	    exit;
	}
	else{
	    $exonblocks->{$gene}={"chr"=>$chr,"strand"=>$strand,"start"=>[],"end"=>[]};

	    my @uniquestart = keys %hashcoords;
	    my @sortedstart = sort {$a <=> $b} @uniquestart;
	    
	    my @start;
	    my @end;
	    my @strand;

	    foreach my $beg (@sortedstart){
		my $en=$hashcoords{$beg};
		push(@start, $beg);
		push(@end, $en);
	    }

	    my $nbex=@start;

	    my $currentstart=$start[0];
	    my $currentend=$end[0];
	 	    
	    for(my $u=1;$u<$nbex;$u++){
		my $thisstart=$start[$u];
		my $thisend=$end[$u];
			
		## cluster blocks if they overlap
		
		if($thisstart>=$currentstart && $thisstart<=($currentend+$collapse)){  
		    
		    ## we only change the end if it's larger than the current position
		    if($thisend>$currentend){
			$currentend=$thisend;
		    }
		}
		else{
		    push(@{$exonblocks->{$gene}{"start"}},$currentstart);
		    push(@{$exonblocks->{$gene}{"end"}},$currentend);
		    
		    $currentstart=$thisstart;
		    $currentend=$thisend;
		 
		}
	    }
	    
	    ## don't forget the last block
	    
	    push(@{$exonblocks->{$gene}{"start"}},$currentstart);
	    push(@{$exonblocks->{$gene}{"end"}},$currentend);
	}
    }
}

##############################################################

sub writeExons{

    my $refblocks=$_[0];
    my $pathoutput=$_[1];

    open(my $output,">".$pathoutput);
    
    foreach my $gene (keys %{$refblocks}){
	my $chr=$refblocks->{$gene}{"chr"};
	my $strand=$refblocks->{$gene}{"strand"};

	my $outstrand=$strand;

	if($strand eq "+"){
	    $outstrand="1";
	}
	
	if($strand eq "-"){
	    $outstrand="-1";
	}
	
	my $nbblocks=@{$refblocks->{$gene}{"start"}};
	
	for(my $i=0;$i<$nbblocks;$i++){
	    my $exonid=$gene.".".($i+1);
	    
	    my $start=${$refblocks->{$gene}{"start"}}[$i];
	    my $end=${$refblocks->{$gene}{"end"}}[$i];
	 
	    print $output $gene."\t".$exonid."\t".$chr."\t".$start."\t".$end."\t".$outstrand."\n";
	}
	
    }
        
    close($output);
}

##############################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script constructs exon blocks from Ensembl annotations. \n";
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
$parameters{"collapseDistance"}=10;
$parameters{"pathOutputExonBlocks"}="NA";

my %defaultvalues;
my @defaultpars=("pathGTF","collapseDistance","pathOutputExonBlocks");


my %numericpars;
my @numericpars=("collapseDistance");

foreach my $par (@numericpars){
    $numericpars{$par}=1;
}

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

#####################################################################
#####################################################################

print "Reading exon coordinates from GTF file...\n";
my %exoncoords;
readGTF($parameters{"pathGTF"},\%exoncoords);
print "Done.\n\n";

my $nbgenes=keys %exoncoords;
print "There are ".$nbgenes." genes.\n";

#####################################################################

print "Making exon blocks...\n";
my %exonblocks;
my $dist=$parameters{"collapseDistance"}+0;

print "Collapsing blocks that are less than ".$dist."bp apart.\n";

makeExonBlocks(\%exoncoords,$dist,\%exonblocks);
print "Done.\n\n";

#####################################################################

print "Writing output...\n";
writeExons(\%exonblocks,$parameters{"pathOutputExonBlocks"});
print "Done.\n\n";

#####################################################################
