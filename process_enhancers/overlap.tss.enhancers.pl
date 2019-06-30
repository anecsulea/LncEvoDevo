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

sub extractTSS{
    my $transcriptcoords=$_[0];
    my $tsscoords=$_[1];
    
    foreach my $tx (keys %{$transcriptcoords}){
	my $chr=$transcriptcoords->{$tx}{"chr"};
	my $strand=$transcriptcoords->{$tx}{"strand"};
	my $start=$transcriptcoords->{$tx}{"start"};
	my $end=$transcriptcoords->{$tx}{"end"};
	my $gene=$transcriptcoords->{$tx}{"gene"};

	if($strand eq "1"){
	    $tsscoords->{$tx}={"gene"=>$gene, "chr"=>$chr, "strand"=>$strand, "position"=>$start};
	} else{
	    if($strand eq "-1"){
		$tsscoords->{$tx}={"gene"=>$gene, "chr"=>$chr, "strand"=>$strand, "position"=>$end};
	    }
	    else{
		print "Weird strand found for tx ".$tx.": ".$strand."\n";
		exit(1);
	    }
	}
    }
}

###########################################################################

sub orderTSS{
    my $coords=$_[0];
    my $ordered=$_[1];
    
    my %hashpos;

    foreach my $tx (keys %{$coords}){
	my $chr=$coords->{$tx}{"chr"};
	my $strand=$coords->{$tx}{"strand"};
	my $position=$coords->{$tx}{"position"};
	
	if(exists $hashpos{$chr}){
	    if(exists $hashpos{$chr}{$position}){
		push(@{$hashpos{$chr}{$position}{"id"}},$tx);
		push(@{$hashpos{$chr}{$position}{"strand"}},$strand);
	    }
	    else{
		$hashpos{$chr}{$position}={"id"=>[$tx], "strand"=>[$strand]};
	    }
	}
	else{
	    $hashpos{$chr}={$position=>{"id"=>[$tx], "strand"=>[$strand]}};
	}
    }

    foreach my $chr (keys %hashpos){
	my @uniqueposition=keys %{$hashpos{$chr}};
	my @sortedposition=sort {$a<=>$b} @uniqueposition;

	foreach my $position (@sortedposition){
	    my $nb=@{$hashpos{$chr}{$position}{"id"}};

	    for(my $i=0; $i<$nb; $i++){
		my $strand=${$hashpos{$chr}{$position}{"strand"}}[$i];

		if(!(exists $ordered->{$chr})){
		    $ordered->{$chr}={"position"=>[], "strand"=>[], "id"=>[]};
		}

		my $id=${$hashpos{$chr}{$position}{"id"}}[$i];
	
		push(@{$ordered->{$chr}{"position"}}, $position);
		push(@{$ordered->{$chr}{"strand"}}, $strand);
		push(@{$ordered->{$chr}{"id"}}, $id);
	    }
	}
    }
}

##############################################################################

sub readEnhancers{
    my $pathin=$_[0];
    my $cgi=$_[1];

    open(my $input, $pathin);

    my $line=<$input>; ## header
    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my $chr=$s[0];
	my $prefix=substr $chr,0,3;

	if($prefix eq "chr"){
	    $chr=substr $chr, 3;
	}

	my $start=$s[1]+0; ## to make it 1-based 
	my $end=$s[2]+0;

	my $id=$chr.",".$start.",".$end;

	$cgi->{$id}={"chr"=>$chr, "start"=>$start, "end"=>$end};
	
	$line=<$input>;
    }

    close($input);
}

##############################################################################

sub orderEnhancersCoords{
    my $coords=$_[0];
    my $ordered=$_[1];
    
    my %hashpos;

    foreach my $id (keys %{$coords}){
	my $chr=$coords->{$id}{"chr"};
	my $start=$coords->{$id}{"start"};
	my $end=$coords->{$id}{"end"};

	if(exists $hashpos{$chr}){
	    if(exists $hashpos{$chr}{$start}){
		push(@{$hashpos{$chr}{$start}{"id"}},$id);
		push(@{$hashpos{$chr}{$start}{"end"}},$end);
	    }
	    else{
		$hashpos{$chr}{$start}={"id"=>[$id],"end"=>[$end]};
	    }
	}
	else{
	    $hashpos{$chr}={$start=>{"id"=>[$id],"end"=>[$end]}};
	}
    }
   

    foreach my $chr (keys %hashpos){
	my @uniquestart=keys %{$hashpos{$chr}};
	my @sortedstart=sort {$a<=>$b} @uniquestart;

	foreach my $start (@sortedstart){
	    my $nb=@{$hashpos{$chr}{$start}{"end"}};

	    for(my $i=0; $i<$nb; $i++){
		my $end=${$hashpos{$chr}{$start}{"end"}}[$i];
		my $id=${$hashpos{$chr}{$start}{"id"}}[$i];
		
		if(!(exists $ordered->{$chr})){
		    $ordered->{$chr}={"start"=>[],"end"=>[], "id"=>[]};
		}

		push(@{$ordered->{$chr}{"start"}}, $start);
		push(@{$ordered->{$chr}{"end"}}, $end);
		push(@{$ordered->{$chr}{"id"}}, $id);
	    }
	}
    }
}

##############################################################################

sub extractOverlapTSSEnhancers{
    my $coords1=$_[0]; ## ordered coordinates
    my $coords2=$_[1]; ## ordered coordinates
    my $margin=$_[2];
    my $overlap=$_[3];

    foreach my $chr (keys %{$coords1}){
	if(exists $coords2->{$chr}){
	    my $nbex1=@{$coords1->{$chr}{"position"}};
	    my $nbex2=@{$coords2->{$chr}{"start"}};
	    
	    my $firstj=0;
	    
	    for(my $i=0; $i<$nbex1; $i++){
		
		my $start1=${$coords1->{$chr}{"position"}}[$i]-$margin;
		my $end1=${$coords1->{$chr}{"position"}}[$i]+$margin;

		my $id1=${$coords1->{$chr}{"id"}}[$i];
			
		my $j=$firstj;
		
		while($j<$nbex2 && ${$coords2->{$chr}{"end"}}[$j]<$start1){
		    $j++;
		}
		
		$firstj=$j;
		
		while($j<$nbex2 && ${$coords2->{$chr}{"start"}}[$j]<=$end1){
		    
		    my $start2=${$coords2->{$chr}{"start"}}[$j];
		    my $end2=${$coords2->{$chr}{"end"}}[$j];
		    my $id2=${$coords2->{$chr}{"id"}}[$j];
		    
		    my $M=max($start1, $start2);
		    my $m=min($end1, $end2);
		    
		    if($M<=$m){
			my $fr1=($m-$M+1)/($end1-$start1+1);
			my $fr2=($m-$M+1)/($end2-$start2+1);
			
			if(exists $overlap->{$id1}){
			    $overlap->{$id1}{$id2}={"start"=>$M,"end"=>$m,"fractionoverlap1"=>$fr1,"fractionoverlap2"=>$fr2};
			}
			else{
			    $overlap->{$id1}={$id2=>{"start"=>$M,"end"=>$m,"fractionoverlap1"=>$fr1,"fractionoverlap2"=>$fr2}};
			}
		    }
		    
		    $j++;
		}
	    }
	}
    }
}


##############################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script computes the overlap between TSS and CpG islands. \n";
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
$parameters{"pathEnhancers"}="NA";
$parameters{"maxDistance"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("pathGTF", "pathEnhancers", "maxDistance", "pathOutput");

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

print "Reading transcript coordinates...\n";

my %exoncoords;
my %transcriptcoords;
my %genetx;

readGTF($parameters{"pathGTF"}, \%exoncoords, \%transcriptcoords, \%genetx);
  
my $nbgene=keys %genetx;
my $nbtx=keys %transcriptcoords;
my $nbex=keys %exoncoords;

print "Found ".$nbgene." genes, ".$nbtx." transcripts and ".$nbex." exons in the annotation.\n";

print "Done.\n";

#####################################################################

print "Extracting and ordering TSS coordinates...\n";

my %tss;
extractTSS(\%transcriptcoords, \%tss);

my %orderedtss;
orderTSS(\%tss, \%orderedtss);

print "Done.\n";

#####################################################################

print "Reading and ordering enhancer coordinates...\n";

my %cgi;

readEnhancers($parameters{"pathEnhancers"}, \%cgi);

my $nbcgi=keys %cgi;

print "There are ".$nbcgi." enhancers.\n";

my %orderedcgi;

orderEnhancersCoords(\%cgi, \%orderedcgi);

print "Done.\n";

#####################################################################

print "Extracting overlap between TSS and enhancers...\n";

my %overlapcgi;

my $margin=$parameters{"maxDistance"}+0;

print "Maximum distance: ".$margin."\n";

extractOverlapTSSEnhancers(\%orderedtss, \%orderedcgi, $margin, \%overlapcgi);
 
print "Done.\n";

#####################################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

print $output "GeneID\tTranscriptID\tChr\tStart\tEnd\tStrand\tTSS\tOverlappingEnhancers\n";

foreach my $gene (keys %genetx){
    foreach my $tx (keys %{$genetx{$gene}}){
	my $chr=$transcriptcoords{$tx}{"chr"};
	my $start=$transcriptcoords{$tx}{"start"};
	my $end=$transcriptcoords{$tx}{"end"};
	my $strand=$transcriptcoords{$tx}{"strand"};
	my $tss=$tss{$tx}{"position"};

	my $idcgi="NA";
	
	if(exists $overlapcgi{$tx}){
	    $idcgi=join(";", keys %{$overlapcgi{$tx}});
	}

	print $output $gene."\t".$tx."\t".$chr."\t".$start."\t".$end."\t".$strand."\t".$tss."\t".$idcgi."\n";
    }
}

close($output);

print "Done.\n";

#####################################################################
