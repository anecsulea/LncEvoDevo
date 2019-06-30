use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

###################################################################
###################################################################

sub readGTF{
    my $pathin=$_[0];
    my $exoncoords=$_[1];
    my $transcriptcoords=$_[2];
    my $genetx=$_[3];
      
    open(my $input, $pathin);
    
    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	
	my $type=$s[2];

	if($type eq "exon"){
	    my $chr=$s[0];
	    my $start=$s[3]+0;
	    my $end=$s[4]+0;
	    my $ss=$s[6];
	    my $strand="NA";
	    
	    if($ss eq "+"){
		$strand="1";
	    } else{
		if($ss eq "-"){
		    $strand="-1";
		} 
	    }

	    ## we discard here unstranded transcripts
	    
	    if($strand ne "NA"){
		my $info=$s[8];
		my @infoarray=split(";", $info);
		
		my $gene=findInfo("gene_id", \@infoarray);
		my $tx=findInfo("transcript_id", \@infoarray);
		
		if($gene eq "NA"){
		    print "Weird! cannot find gene id for ".$line."\n";
		    exit(1);
		}
		
		if($tx eq "NA"){
		    print "Weird! cannot find transcript id for ".$line."\n";
		    exit(1);
		}
		
		if(exists $genetx->{$gene}){
		    $genetx->{$gene}{$tx}=1;
		} else{
		    $genetx->{$gene}={$tx=>1};
		}

		my $exonid=$chr.",".$start.",".$end.",".$strand;
		
		$exoncoords->{$exonid}={"start"=>$start, "end"=>$end, "chr"=>$chr, "strand"=>$strand};
		
		if(exists $transcriptcoords->{$tx}){
		    if($start < $transcriptcoords->{$tx}{"start"}){
			$transcriptcoords->{$tx}{"start"}=$start;
		    }
		    
		    if($end > $transcriptcoords->{$tx}{"end"}){
			$transcriptcoords->{$tx}{"end"}=$end;
		    }
		    
		} else{
		    $transcriptcoords->{$tx}={"start"=>$start, "end"=>$end, "chr"=>$chr, "strand"=>$strand, "gene"=>$gene};
		}
	    }
	}
	
	$line=<$input>;
    }
    
    close($input);
}

#########################################################################

sub findInfo{
    my $pattern=$_[0];
    my $array=$_[1];
    
    my $res="NA";
    my $nbfound=0;
    
    my @grepres=grep(/${pattern}/,@{$array});
    
    if(@grepres==1){
	my @t=split("\"",$grepres[0]);
	$res=$t[1];
	$nbfound++;
    }
    else{
	foreach my $g (@grepres){
	    my @u=split(" ",$g);
	    
	    if($u[0] eq $pattern){
		$nbfound++;
		my @t=split("\"",$g);
		$res=$t[1];
	    }
	}
    }
    
    if($nbfound==1){
	return $res;
    } else{
	return "NA";
    }
}

###########################################################################

sub readFasta{
    my $path=$_[0];
    my $reffasta=$_[1];
   
    my @s=split("\\.",$path);
    my $ext=$s[-1];

    my $input;

    if($ext eq "gz"){
	open($input,"zcat $path |");
    }
    else{
	open($input, $path);
    }
    
    my $line=<$input>;

    while($line){
	my $b=substr $line,0,1;
	
	if($b eq ">"){
	    chomp $line;
	    my $id=substr $line,1;

	    my @s=split(" ",$id);
	    $id=$s[0];

	    # print "saw chromosome ".$id."\n";
	    
	    $reffasta->{$id}="";

	    $line=<$input>;
	    $b=substr $line,0,1;
	    
	    while($line && !($b eq ">")){
		chomp $line;
		$reffasta->{$id}.=$line;
		$line=<$input>;
		$b=substr $line,0,1;
	    }
	}
    }

    close($input);
}

############################################################################

sub reverseComplement{
    my $sequence=$_[0];
    
    my $rev=reverse $sequence;

    $rev=~s/A/X/g;
    $rev=~s/C/Y/g;
    $rev=~s/G/Z/g;
    $rev=~s/T/W/g;

    $rev=~s/X/T/g;
    $rev=~s/Y/G/g;
    $rev=~s/Z/C/g;
    $rev=~s/W/A/g;

    return $rev;
}


###############################################################

sub writeSequence{
    my $sequence=$_[0];
    my $name=$_[1];
    my $output=$_[2];

    my $n=length $sequence;

    print $output ">".$name."\n";

    my $i=0;

    while($i<($n-60)){

        my $subseq=substr $sequence,$i,60;

        print $output $subseq ."\n";

        $i+=60;
    }

    if($i<$n){
        my $subseq=substr $sequence,$i;
        print $output $subseq ."\n";
    }
}

##############################################################

sub readChromosomeCorrespondence{
    my $pathin=$_[0]; 
    my $corresp12=$_[1]; 
    my $corresp21=$_[2]; 

    open(my $input, $pathin);

    my $line=<$input>; ## header
    
    my %dupliens;
    my %dupliucsc;
    
    while($line){
	chomp $line;
	my @s=split("\t", $line);
	
	my $ens=$s[0];
	my $ucsc=$s[1];
	
	if(exists $corresp12->{$ens}){
	    $dupliens{$ens}=1;
	    $dupliucsc{$ucsc}=1;
	}
	else{
	    $corresp12->{$ens}=$ucsc;
	}
	
	if(exists $corresp21->{$ucsc}){
	    $dupliens{$ens}=1;
	    $dupliucsc{$ucsc}=1;
	}
	else{
	    $corresp21->{$ucsc}=$ens;
	}
	
	$line=<$input>;
    }
    
    close($input);
    
    my $nbdupliens=keys %dupliens;
    
    if($nbdupliens>0){
	print "Found ".$nbdupliens." ambiguous Ensembl chromosomes: ".join("; ",keys %dupliens).".\n";

	foreach my $idens (keys %dupliens){
	    delete $corresp12->{$idens};
	}
    }

    my $nbdupliucsc=keys %dupliucsc;
    
    if($nbdupliucsc>0){
	print "Found ".$nbdupliucsc." ambiguous UCSC chromosomes: ".join("; ", keys %dupliucsc).".\n";

	foreach my $iducsc (keys %dupliucsc){
	    delete $corresp21->{$iducsc};
	}
    }
}

###########################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script extracts promoter sequences and coordinates. \n";
    print "\n";
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

############################################################################
############################################################################

my %parameters;

$parameters{"pathGTF"}="NA";
$parameters{"pathGenomeSequence"}="NA";
$parameters{"promoterSize"}="NA";
$parameters{"pathChromosomeCorrespondence"}="NA";
$parameters{"pathOutputCoords"}="NA";
$parameters{"pathOutputSequences"}="NA";

my @defaultpars=("pathGTF", "pathGenomeSequence", "promoterSize", "pathChromosomeCorrespondence", "pathOutputCoords", "pathOutputSequences");

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

##############################################################################
##############################################################################

print "Reading annotations...\n";

my %exoncoords;
my %transcriptcoords;
my %genetx;

readGTF($parameters{"pathGTF"}, \%exoncoords, \%transcriptcoords, \%genetx);

my $nbex=keys %exoncoords;
my $nbtx=keys %transcriptcoords;
my $nbg=keys %genetx;

print "Found ".$nbex." exons, ".$nbtx." transcripts and ".$nbg." genes.\n";
print "Done.\n";

##############################################################################

print "Reading genome sequence...\n";

my %genome;
readFasta($parameters{"pathGenomeSequence"}, \%genome);

my $nbchr=keys %genome;

print "Saw ".$nbchr." chromosomes.\n";

print "Done.\n";

##############################################################################

print "Reading chromosome correspondence...\n";

my %ensucsc;
my %ucscens;

readChromosomeCorrespondence($parameters{"pathChromosomeCorrespondence"}, \%ensucsc, \%ucscens);

print "Done.\n";

##############################################################################

print "Writing output for promoter coordinates...\n";

open(my $outputfasta, ">".$parameters{"pathOutputSequences"});
open(my $output, ">".$parameters{"pathOutputCoords"});

my $promsize=$parameters{"promoterSize"}+0;

print "Promoter size: ".$promsize."bp.\n"; 

foreach my $gene (keys %genetx){
    my %tss;
    my $chr;
    my $strand;
    
    foreach my $tx (keys %{$genetx{$gene}}){
	$chr=$transcriptcoords{$tx}{"chr"};
	$strand=$transcriptcoords{$tx}{"strand"};

	my $start=$transcriptcoords{$tx}{"start"};
	my $end=$transcriptcoords{$tx}{"end"};
       
	my $tsscoords=$chr.",".($start-$promsize).",".($start-1).",".$strand;

	if($strand eq "-1"){
	    $tsscoords=$chr.",".($end+1).",".($end+$promsize).",".$strand;
	}

	if(exists $tss{$tsscoords}){
	    $tss{$tsscoords}{$tx}=1;
	} else{
	    $tss{$tsscoords}={$tx=>1};
	}
    }

    foreach my $coord (keys %tss){
	my $id=$gene."_".join(";", keys %{$tss{$coord}});
	my @s=split(",", $coord);
	my $promstart=$s[1];
	my $promend=$s[2];
	
	my $newstrand="NA";

	if($strand eq "1"){
	    $newstrand="+";
	}
	else{
	    if($strand eq "-1"){
		$newstrand="-";
	    } else{
		print "Weird strand for ".$gene."\n";
		exit(1);
	    }
	}

	if(exists $ensucsc{$chr}){
	    my $newchr=$ensucsc{$chr};
	    
	    print $output $newchr."\t".($promstart-1)."\t".$promend."\t".$id."\t1000\t".$newstrand."\n";
	    
	} else{
	    print "Cannot find ".$chr." in Ensembl-UCSC correspondence.\n";
	}

	if(exists $genome{$chr}){
	    my $seq=substr $genome{$chr}, ($promstart-1), ($promend-$promstart+1);

	    if($strand eq "-1"){
		$seq=reverseComplement($seq);
	    }

	    writeSequence($seq, $id, $outputfasta);
	}
	else{
	    print "Weird! cannot find chromosome ".$chr." in genome sequence.\n";
	    exit(1);
	}
    }
}

close($output);
close($outputfasta);

print "Done.\n";

##############################################################################

