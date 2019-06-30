use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

###################################################################
###################################################################

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

###################################################################

sub readGeneInfo{
    my $pathInfo=$_[0];
    my $info=$_[1];
    
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
	my $thisbio=$s[$header{"biotype"}];
	my $descr=$s[$header{"description"}];

	$info->{$id}={"biotype"=>$thisbio,"description"=>$descr};

	$line=<$input>;
    }
    
    close($input);
}

###################################################################

sub readSynonyms{
    my $pathin=$_[0];
    my $oldnew=$_[1];
    my $newold=$_[2];

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

	my $old=$s[$header{"OldGeneID"}];
	my $new=$s[$header{"NewGeneID"}];
	
	my @oldgenes=split(",", $old);
	
	foreach my $og (@oldgenes){
	    if(exists $oldnew->{$og}){
		$oldnew->{$og}{$new}=1;
	    }
	    else{
		$oldnew->{$og}={$new=>1};
	    }
 
	    if(exists $newold->{$new}){
		$newold->{$new}{$og}=1;	
	    }
	    else{
		$newold->{$new}={$og=>1};	
	    }
	}

	$line=<$input>;
    }

    close($input);
}

###################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script extracts intron block coordinates. \n";
    print "\n";
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

######################################################################
######################################################################

my %parameters;

$parameters{"pathExonBlocks"}="NA";
$parameters{"pathGeneInfo"}="NA";
$parameters{"pathSynonyms"}="NA";
$parameters{"biotypes"}="NA";
$parameters{"pathOutputIntronBlocks"}="NA";
$parameters{"pathOutputSelectedExonBlocks"}="NA";

my %defaultvalues;
my @defaultpars=("pathExonBlocks", "pathGeneInfo", "pathSynonyms", "biotypes", "pathOutputIntronBlocks", "pathOutputSelectedExonBlocks");

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

#####################################################################
#####################################################################

print "Reading exon blocks...\n";

my %exonblocks;

readExonBlocks($parameters{"pathExonBlocks"}, \%exonblocks);

my $nballgenes=keys %exonblocks;

print "There are ".$nballgenes." genes.\n";

print "Done.\n";

#####################################################################

print "Reading gene info...\n";

my %geneinfo;

readGeneInfo($parameters{"pathGeneInfo"}, \%geneinfo);

print "Done.\n";

#####################################################################

print "Reading synonyms...\n";

my %oldnew;
my %newold;

readSynonyms($parameters{"pathSynonyms"}, \%oldnew, \%newold);

my $nbn=keys %newold;

print $nbn." new genes in synonyms.\n";

print "Done.\n";

#####################################################################

print "Selecting genes depending on their biotypes...\n";

my @bio=split(",", $parameters{"biotypes"});

my %accbio;

foreach my $b (@bio){
    $accbio{$b}=1;
}

my $nb=@bio;

print "There are ".$nb." accepted biotypes: ".join(", ", @bio)."\n";

my %selectedgenes;

my $nborig=0;
my $nbsyn=0;

foreach my $gene (keys %exonblocks){
    if(exists $accbio{"all"}){
	$selectedgenes{$gene}=1;
	$nborig++;
    }
    else{
	if(exists $geneinfo{$gene}){
	    my $thisbio=$geneinfo{$gene}{"biotype"};
	    if(exists $accbio{$thisbio}){
		$selectedgenes{$gene}=1;
		$nborig++;
	    }
	} else{
	    if(exists $newold{$gene}){
		foreach my $oldid (keys %{$newold{$gene}}){
		    if(exists $geneinfo{$oldid}){
			my $thisbio=$geneinfo{$oldid}{"biotype"};
			
			if(exists $accbio{$thisbio}){
			    $selectedgenes{$gene}=1;
			    $nbsyn++;
			    last;
			}  
		    }
		}
	    }
	} 
    }
}

my $nbsel=keys %selectedgenes;

print "We selected ".$nbsel." genes with the correct biotypes: ".$nborig." in the original annotations, ".$nbsyn." based on synonyms.\n";

print "Done.\n";

#####################################################################

print "Extracting intron blocks and writing output...\n";

open(my $output, ">".$parameters{"pathOutputIntronBlocks"});

foreach my $gene (keys %selectedgenes){
    my $nbex=@{$exonblocks{$gene}{"start"}};

    if($nbex>=2){
	my $chr=$exonblocks{$gene}{"chr"};
	my $strand=$exonblocks{$gene}{"strand"};

	for(my $i=0; $i<($nbex-1); $i++){
	    my $endex=${$exonblocks{$gene}{"end"}}[$i];
	    my $startex=${$exonblocks{$gene}{"start"}}[$i+1];

	    if(($endex+1) <= ($startex-1)){
		print $output $gene."\t".$gene.".".($i+1)."\t".$chr."\t".($endex+1)."\t".($startex-1)."\t".$strand."\n";
	    } else{
		if($endex > $startex){
		    print "Weird! exon blocks not ordered for ".$gene."\n";
		    print $startex." ".$endex."\n";
		    exit(1);
		}
	    }
	}
    }
}

close($output);

print "Done.\n";


#####################################################################

if($nbsel < $nballgenes){
    print "Writing output also for selected exon blocks.\n";
    
    open(my $output, ">".$parameters{"pathOutputSelectedExonBlocks"});

    foreach my $gene (keys %selectedgenes){
	my $nbex=@{$exonblocks{$gene}{"start"}};
	my $chr=$exonblocks{$gene}{"chr"};
	my $strand=$exonblocks{$gene}{"strand"};

	for(my $i=0; $i<$nbex; $i++){
	    my $startex=${$exonblocks{$gene}{"start"}}[$i];
	    my $endex=${$exonblocks{$gene}{"end"}}[$i];

	    print $output $gene."\t".$gene.".".($i+1)."\t".$chr."\t".$startex."\t".$endex."\t".$strand."\n";
	}
    }

    close($output);

    print "Done.\n";
}

#####################################################################
