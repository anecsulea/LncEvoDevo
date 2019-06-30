use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################

sub readExonCoords{

    my $pathExons=$_[0];
    my $exons=$_[1];

    open(my $input,$pathExons);
    
    my $line=<$input>; ## header
    my %header;
    chomp $line;
    my @s=split("\t",$line);
    
    for(my $i=0;$i<@s;$i++){
	$header{$s[$i]}=$i;
    }

    $line=<$input>;
    
    while($line){
	chomp $line;
	my @s=split("\t",$line);
	
	my $id=$s[$header{"stable_id"}];
	my $chr=$s[$header{"name"}];
	my $start=$s[$header{"seq_region_start"}];
	my $end=$s[$header{"seq_region_end"}];
	my $strand=$s[$header{"seq_region_strand"}];
	
	$exons->{$id}={"chr"=>$chr,"start"=>$start,"end"=>$end,"strand"=>$strand};

	$line=<$input>;
    }

    close($input);
}


##############################################################

sub extractReadThroughGenes{
    my $pathInfo=$_[0];
    my $readthrough=$_[1];

    open(my $input,$pathInfo);
    
    my $line=<$input>; ## header
    my %header;
    chomp $line;
    my @s=split("\t",$line);
    
    for(my $i=0;$i<@s;$i++){
	$header{$s[$i]}=$i;
    }
    
    $line=<$input>;
    
    while($line){
	chomp $line;
	my @s=split("\t",$line);
	
	my $id=$s[$header{"stable_id"}];

	$readthrough->{$id}=1;
	
	$line=<$input>;
    }
    
    close($input);

}


##############################################################

sub extractReadThroughTranscripts{
    my $pathInfo=$_[0];
    my $readthrough=$_[1];

    open(my $input,$pathInfo);
    
    my $line=<$input>; ## header
    my %header;
    chomp $line;
    my @s=split("\t",$line);
    
    for(my $i=0;$i<@s;$i++){
	$header{$s[$i]}=$i;
    }
    
    $line=<$input>;
    
    while($line){
	chomp $line;
	my @s=split("\t",$line);
	
	my $id=$s[$header{"TranscriptID"}];

	$readthrough->{$id}=1;
	
	$line=<$input>;
    }
    
    close($input);

}

##############################################################

sub readExonAssignment{

    my $pathin=$_[0];
    my $txex=$_[1];
  
    open(my $input,$pathin);
    
    my $line=<$input>; ## header
    $line=<$input>;

    while($line){
	chomp $line;
	
	my @s=split("\t",$line);
	my $exid=$s[0];
	my $txid=$s[1];

	if(exists $txex->{$txid}){
	    push(@{$txex->{$txid}},$exid);
	}
	else{
	    $txex->{$txid}=[$exid];
	}

	
	$line=<$input>;
    }

    close($input);
}

##############################################################

sub readGeneBiotype{

    my $pathBiotype=$_[0];
    my $bio=$_[1];
    
    open(my $input,$pathBiotype);
    
    my $line=<$input>; ## header
    my %header;
    chomp $line;
    my @s=split("\t",$line);
    
    for(my $i=0;$i<@s;$i++){
	$header{$s[$i]}=$i;
    }
    
    $line=<$input>;
    
    while($line){
	chomp $line;
	my @s=split("\t",$line);
	
	my $id=$s[$header{"stable_id"}];
	my $thisbio=$s[$header{"biotype"}];

	$bio->{$id}=$thisbio;

	$line=<$input>;
    }
    
    close($input);
    
}

##############################################################

sub readTranscriptInfo{

    my $pathBiotype=$_[0];
    my $bio=$_[1];
    my $genetx=$_[2];

    open(my $input,$pathBiotype);
    
    my $line=<$input>; ## header
    my %header;
    chomp $line;
    my @s=split("\t",$line);
    
    for(my $i=0;$i<@s;$i++){
	$header{$s[$i]}=$i;
    }

    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t",$line);
	
	my $geneid=$s[0];
	my $id=$s[1];
	my $thisbio=$s[$header{"biotype"}];

	$bio->{$id}=$thisbio;

	if(exists $genetx->{$geneid}){
	    push(@{$genetx->{$geneid}},$id);
	}
	else{
	    $genetx->{$geneid}=[$id];
	}
	
	$line=<$input>;
    }
    
    close($input);
}

##############################################################

sub makeExonBlocks{

    my $genetx=$_[0];
    my $txex=$_[1];
    my $genebio=$_[2];
    my $txbio=$_[3];
    my $excoord=$_[4];
    my $collapse=$_[5];
    my $filter=$_[6];
    my $readthrough=$_[7];
    my $exonblocks=$_[8];

    foreach my $gene (keys %{$genetx}){
	my $thisbiogene=$genebio->{$gene};
	
	my %hashbegin;
	my $chr="NA";

	foreach my $tx (@{$genetx->{$gene}}){
	    my $thisbiotx=$txbio->{$tx};

	    if(($filter eq "yes" && ($thisbiogene eq "protein_coding" && $thisbiotx eq "protein_coding") || ($thisbiogene ne "protein_coding")) || ($filter eq "no")){
		
		if(!(exists $readthrough->{$tx})){
		    foreach my $ex (@{$txex->{$tx}}){
			if($chr eq "NA"){
			    $chr=$excoord->{$ex}{"chr"};
			}
			
			my $strand=$excoord->{$ex}{"strand"};
			my $b=$excoord->{$ex}{"start"};
			my $e=$excoord->{$ex}{"end"};
			
			if(exists $hashbegin{$b}){
			    if(exists $hashbegin{$b}{$e}){
				push(@{$hashbegin{$b}{$e}}, $strand);
			    }
			    else{
				$hashbegin{$b}{$e}=[$strand];
			    }
			}
			else{
			    $hashbegin{$b}={$e=>[$strand]};
			}
		    }
		}
	    }
	}

	if($chr eq "NA"){
	    print $gene." doesn't have any filtered transcripts\n";
	}
	else{
	    $exonblocks->{$gene}={"chr"=>$chr,"strand"=>[],"begin"=>[],"end"=>[]};

	    my @uniquebegin = keys %hashbegin;
	    my @sortedbegin = sort {$a <=> $b} @uniquebegin;
	    
	    my @begin;
	    my @end;
	    my @strand;

	    foreach my $beg (@sortedbegin){
		my @thisend=keys %{$hashbegin{$beg}};
		my @sortedend = sort {$a <=> $b} @thisend;

		foreach my $en (@sortedend){
		    foreach my $str (@{$hashbegin{$beg}{$en}}){
			push(@begin, $beg);
			push(@end, $en);
			push(@strand, $str);
		    }
		}
	    }

	    my $nbex=@begin;

	    my $currentbegin=$begin[0];
	    my $currentend=$end[0];
	    my $currentstrand=$strand[0];
	    
	    for(my $u=1;$u<$nbex;$u++){
		my $thisbegin=$begin[$u];
		my $thisend=$end[$u];
		my $thisstrand=$strand[$u];
		
		## cluster blocks if they overlap
		
		if($thisbegin>=$currentbegin && $thisbegin<=($currentend+$collapse) && $thisstrand eq $currentstrand){  
		    
		    ## we only change the end if it's larger than the current position
		    if($thisend>$currentend){
			$currentend=$thisend;
		    }
		}
		else{
		    push(@{$exonblocks->{$gene}{"begin"}},$currentbegin);
		    push(@{$exonblocks->{$gene}{"end"}},$currentend);
		    push(@{$exonblocks->{$gene}{"strand"}},$currentstrand);
		    
		    $currentbegin=$thisbegin;
		    $currentend=$thisend;
		    $currentstrand=$thisstrand;
		}
		
	    }
	    
	    ## don't forget the last block
	    
	    push(@{$exonblocks->{$gene}{"begin"}},$currentbegin);
	    push(@{$exonblocks->{$gene}{"end"}},$currentend);
	    push(@{$exonblocks->{$gene}{"strand"}},$currentstrand);
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
	
	my $nbblocks=@{$refblocks->{$gene}{"begin"}};
	
	for(my $i=0;$i<$nbblocks;$i++){
	    my $exonid=$gene.".".($i+1);
	    
	    my $begin=${$refblocks->{$gene}{"begin"}}[$i];
	    my $end=${$refblocks->{$gene}{"end"}}[$i];
	    my $strand=${$refblocks->{$gene}{"strand"}}[$i];
	    
	    print $output $gene."\t".$exonid."\t".$chr."\t".$begin."\t".$end."\t".$strand."\n";
	}
	
    }
        
    close($output);
}

##############################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script constructs exon blocks from Ensembl annotations. Transcripts annotated as retained_intron are excluded.\n";
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

$parameters{"pathExonCoords"}="NA";
$parameters{"pathExonAssignment"}="NA";
$parameters{"pathGeneInfo"}="NA";
$parameters{"pathReadThrough"}="NA";
$parameters{"pathTranscriptInfo"}="NA";
$parameters{"collapseDistance"}=10;
$parameters{"filter"}="NA";
$parameters{"pathOutputExonBlocks"}="NA";

my %defaultvalues;
my @defaultpars=("pathExonCoords","pathExonAssignment","pathGeneInfo","pathTranscriptInfo","pathReadThrough","collapseDistance","filter","pathOutputExonBlocks");


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

print "Reading exon coordinates...\n";
my %exoncoords;
readExonCoords($parameters{"pathExonCoords"},\%exoncoords);
print "Done.\n\n";

my $nbexons=keys %exoncoords;
print "There are ".$nbexons." exons.\n";

print "Reading exon-transcript assignment...\n";
my %transcriptexons;
readExonAssignment($parameters{"pathExonAssignment"},\%transcriptexons);
print "Done.\n\n";

print "Reading gene biotypes...\n";
my %genebio;
readGeneBiotype($parameters{"pathGeneInfo"},\%genebio);
print "Done.\n\n";

print "Reading transcript info...\n";
my %txbio;
my %genetranscripts;
readTranscriptInfo($parameters{"pathTranscriptInfo"},\%txbio,\%genetranscripts);
print "Done.\n\n";

my $nbgenes=keys %genetranscripts;
print "There are ".$nbgenes." genes.\n";

my $nbtx=keys %txbio;
print "There are ".$nbtx." transcripts.\n";

#####################################################################

print "Identifying read-through genes...\n";
my %readthrough;
extractReadThroughTranscripts($parameters{"pathReadThrough"},\%readthrough);
my $nbr=keys %readthrough;

print "There are ".$nbr." readthrough transcripts.\n";
print "Done.\n\n";

#####################################################################

print "Making exon blocks...\n";
my %exonblocks;
print "Collapsing blocks that are less than ".$parameters{"collapseDistance"}."bp apart.\n";
my $filter=$parameters{"filter"};
if($filter ne "yes" && $filter ne "no"){
    print "unknown filter!!! ".$filter."\n";
}
print "filter ".$filter."\n";
makeExonBlocks(\%genetranscripts,\%transcriptexons,\%genebio, \%txbio,\%exoncoords,$parameters{"collapseDistance"},$filter, \%readthrough, \%exonblocks);
print "Done.\n\n";

#####################################################################

print "Removing readthrough genes...\n";

foreach my $gene (keys %readthrough){
    delete $exonblocks{$gene};
}

my $nbg=keys %exonblocks;

print $nbg." genes remaining.\n";
print "Done.\n";

#####################################################################

print "Writing output...\n";
writeExons(\%exonblocks,$parameters{"pathOutputExonBlocks"});
print "Done.\n\n";

#####################################################################
