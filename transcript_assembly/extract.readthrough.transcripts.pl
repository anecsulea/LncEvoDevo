#!/usr/bin/perl

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################
sub  trim {
    my $s = $_[0]; 
    $s =~ s/^\s+|\s+$//g;
    return $s
}

##############################################################

sub readGTF{
    my $pathin=$_[0];
    my $exoncoords=$_[1];
    my $exontx=$_[2];
    my $exongene=$_[3];
    my $genetx=$_[4];
    my $txex=$_[5];
    
    open(my $input, $pathin);
    
    my $line=<$input>;
    my $prefix=substr $line, 0,1;
    
    while($prefix eq "#"){
	$line=<$input>;
	$prefix=substr $line, 0,1;
    }

    my $nbdiscarded=0; ## we discard genes with "." strand

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


		my $geneid=findInfo("^gene_id", \@t);
		my $txid=findInfo(" transcript_id", \@t);
		
	
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
		
		$exoncoords->{$exonid}={"chr"=>$chr,"start"=>$start, "end"=>$end, "strand"=>$strand};
		
		## exon - tx correspondence
		
		if(exists $exontx->{$exonid}){
		    $exontx->{$exonid}{$txid}=1;
		}
		else{
		    $exontx->{$exonid}={$txid=>1};
		}
		
		if(exists $exongene->{$exonid}){
		    $exongene->{$exonid}{$geneid}=1;
		}
		else{
		    $exongene->{$exonid}={$geneid=>1};
		}
		
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
    
    my @grepres=grep(/^${pattern}/,@{$array});

    my $nbg=@grepres;

    if(@grepres==1){
	my @t=split("\"",$grepres[0]);
	$res=$t[1];
    }
    
    return $res;
}

##############################################################

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

##############################################################

sub writeExonBlocks{
    my $refexons=$_[0];
    my $synonyms=$_[1];
    my $pathout=$_[2];
    
    open(my $output,">".$pathout);

    foreach my $gene (keys %{$refexons}){
	my $chr=$refexons->{$gene}{"chr"};
	my $strand=$refexons->{$gene}{"strand"};
	my $nbexons=@{$refexons->{$gene}{"start"}};

	my $syn="NA";

	if(exists $synonyms->{$gene}){
	    $syn=$synonyms->{$gene};
	}
	
	for(my $i=0;$i<$nbexons;$i++){
	    my $start=${$refexons->{$gene}{"start"}}[$i];
	    my $end=${$refexons->{$gene}{"end"}}[$i];
	    
	    my $idexon=$gene.".".($i+1);
	    print $output $gene."\t".$idexon."\t".$chr."\t".$start."\t".$end."\t".$strand."\t".$syn."\n";
	}
    }
    
    close($output);
}


##############################################################

sub orderExons{
    my $exons=$_[0];
    my $ordered=$_[1];
    
    my %hashpos;

    foreach my $exid (keys %{$exons}){
	my $chr=$exons->{$exid}{"chr"};
	my $strand=$exons->{$exid}{"strand"};
	my $start=$exons->{$exid}{"start"};
	my $end=$exons->{$exid}{"end"};

	if(exists $hashpos{$chr}){
	    if(exists $hashpos{$chr}{$start}){
		push(@{$hashpos{$chr}{$start}{"id"}},$exid);
		push(@{$hashpos{$chr}{$start}{"end"}},$end);
		push(@{$hashpos{$chr}{$start}{"strand"}},$strand);
	    }
	    else{
		$hashpos{$chr}{$start}={"id"=>[$exid],"end"=>[$end],"strand"=>[$strand]};
	    }
	}
	else{
	    $hashpos{$chr}={$start=>{"id"=>[$exid],"end"=>[$end],"strand"=>[$strand]}};
	}
    }

    foreach my $chr (keys %hashpos){
	$ordered->{$chr}={"start"=>[],"end"=>[],"strand"=>[],"id"=>[]};

	my @uniquestart=keys %{$hashpos{$chr}};
	my @sortedstart=sort {$a<=>$b} @uniquestart;

	foreach my $start (@sortedstart){
	    my $nb=@{$hashpos{$chr}{$start}{"end"}};

	    for(my $i=0; $i<$nb; $i++){
		my $end=${$hashpos{$chr}{$start}{"end"}}[$i];
		my $strand=${$hashpos{$chr}{$start}{"strand"}}[$i];
		my $id=${$hashpos{$chr}{$start}{"id"}}[$i];

		push(@{$ordered->{$chr}{"start"}}, $start);
		push(@{$ordered->{$chr}{"end"}}, $end);
		push(@{$ordered->{$chr}{"strand"}}, $strand);
		push(@{$ordered->{$chr}{"id"}}, $id);
	    }
	}
    }
}

##############################################################

sub extractOverlap{
    my $coords1=$_[0]; ## ordered coordinates
    my $coords2=$_[1]; ## ordered coordinates
    my $margin=$_[2];
    my $type=$_[3];
    my $overlap=$_[4];

    foreach my $chr (keys %{$coords1}){
	if(exists $coords2->{$chr}){
	    my $nbex1=@{$coords1->{$chr}{"start"}};
	    my $nbex2=@{$coords2->{$chr}{"start"}};
	    
	    my $firstj=0;
	    
	    for(my $i=0; $i<$nbex1; $i++){
		
		my $start1=${$coords1->{$chr}{"start"}}[$i]-$margin;
		my $end1=${$coords1->{$chr}{"end"}}[$i]+$margin;
		my $strand1=${$coords1->{$chr}{"strand"}}[$i];
		my $id1=${$coords1->{$chr}{"id"}}[$i];
		
		my $j=$firstj;
		
		while($j<$nbex2 && ${$coords2->{$chr}{"end"}}[$j]<$start1){
		    $j++;
		}
		
		$firstj=$j;
		
		while($j<$nbex2 && ${$coords2->{$chr}{"start"}}[$j]<=$end1){
		    
		    my $start2=${$coords2->{$chr}{"start"}}[$j];
		    my $end2=${$coords2->{$chr}{"end"}}[$j];
		    my $strand2=${$coords2->{$chr}{"strand"}}[$j];
		    my $id2=${$coords2->{$chr}{"id"}}[$j];
		    
		    if(($strand1 eq $strand2 && $type eq "sense") || ($strand1 ne $strand2 && $type eq "antisense")){
			
			my $M=max($start1,$start2);
			my $m=min($end1,$end2);
			
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
		    }
		    
		    $j++;
		}
	    }
	}
    }
}


##########################################################################

sub addToCluster{
    my $refconnected=$_[0];
    my $refclusters=$_[1];
    my $indexclusters=$_[2];
    my $key=$_[3];

    if(exists $refconnected->{$key}){
                    
        my $indexcluster="NA";
	
        if(exists $indexclusters->{$key}){
            $indexcluster=$indexclusters->{$key}
        }       
              
        if($indexcluster eq "NA"){
            my $nbclusters=keys %{$refclusters};
            $indexcluster=$nbclusters+1;
            $refclusters->{$indexcluster}=[$key];
            $indexclusters->{$key}=$indexcluster;
        }
                
        foreach my $connection (keys %{$refconnected->{$key}}){
            if(!(exists $indexclusters->{$connection})){
                push(@{$refclusters->{$indexcluster}},$connection);
                $indexclusters->{$connection}=$indexcluster;
                addToCluster($refconnected,$refclusters,$indexclusters,$connection);
            }
        }

        delete $refconnected->{$key};
    }
}

###########################################################################

sub extractClusters{
    my $refconnected=$_[0];
    my $refclusters=$_[1];
    
    my %indexclusters;
    
    my $nbconnected=keys %{$refconnected};
    
    my $round=0;
    
    my %alreadyprinted;
    
    while($nbconnected>0){
	
        foreach my $key (keys %{$refconnected}){
            addToCluster($refconnected,$refclusters,\%indexclusters,$key);
            $nbconnected=keys %{$refconnected};
        }
        
        $round++;
        
        $nbconnected=keys %{$refconnected};
    }
}

##########################################################################

sub makeExonBlocks{
    my $genes=$_[0];
    my $exons=$_[1];
    my $collapse=$_[2];
    my $exonblocks=$_[3];

    foreach my $gene (keys %{$genes}){
	my $chr="NA";
	my $strand="NA";
	
	my %hashbegin;
		
	foreach my $exon (@{$genes->{$gene}}){
	    if($chr eq "NA"){
		$chr=$exons->{$exon}{"chr"};
		$strand=$exons->{$exon}{"strand"};
	    }
	    else{
		if($exons->{$exon}{"chr"} ne $chr || $exons->{$exon}{"strand"} ne $strand){
		    print "Weird!!! saw chromosomes ".$exons->{$exon}{"chr"}." and ". $chr.", strands ". $exons->{$exon}{"strand"}." and ". $strand." for ".$gene."\n";
		    exit(1);
		}
	    }
	    
	    my $b=$exons->{$exon}{"start"};
	    my $e=$exons->{$exon}{"end"};
	    
	    if(exists $hashbegin{$b}){
		$hashbegin{$b}{$e}=1;
	    }
	    else{
		$hashbegin{$b}={$e=>1};
	    }
	}	    
	  

	$exonblocks->{$gene}={"chr"=>$chr,"strand"=>$strand,"start"=>[],"end"=>[]};
	
	### now make exon blocks

	my @uniquebegin = keys %hashbegin;
	my @sortedbegin = sort {$a <=> $b} @uniquebegin;
	
	my @begin;
	my @end;
	
	foreach my $beg (@sortedbegin){
	    my @thisend=keys %{$hashbegin{$beg}};
	    my @sortedend = sort {$a <=> $b} @thisend;
	    
	    foreach my $en (@sortedend){
		push(@begin, $beg);
		push(@end, $en);
	    }
	}	
	
	my $nbex=@begin;
	
	my $currentbegin=$begin[0];
	my $currentend=$end[0];
	
	for(my $u=1;$u<$nbex;$u++){
	    my $thisbegin=$begin[$u];
	    my $thisend=$end[$u];
	    		
	    ## cluster blocks if they overlap
	    
	    if($thisbegin>=$currentbegin && $thisbegin<=($currentend+$collapse)){  
		
		## we only change the end if it's larger than the current position
		if($thisend>$currentend){
		    $currentend=$thisend;
		}
	    }
	    else{
		push(@{$exonblocks->{$gene}{"start"}},$currentbegin);
		push(@{$exonblocks->{$gene}{"end"}},$currentend);
		
		$currentbegin=$thisbegin;
		$currentend=$thisend;
	    }
	}
	
	## don't forget the last block
	
	push(@{$exonblocks->{$gene}{"start"}},$currentbegin);
	push(@{$exonblocks->{$gene}{"end"}},$currentend);
    }
}

##############################################################

sub extractChimericTranscripts{
    my $txexons=$_[0];
    my $exongene=$_[1];
    my $monogenes=$_[2];
    my $exonoverlap=$_[3];
    my $chimerictx=$_[4];
  
    foreach my $tx (keys %{$txexons}){	
	my $prefix=substr $tx, 0, 3; ## only assembled transcripts

	my %ensgenes; ## Ensembl genes with which this transcript overlaps (including its own)
	
	foreach my $exon (keys %{$txexons->{$tx}}){
	    
	    my %ensgenesex;
	    
	    if(exists $exonoverlap->{$exon}){
		foreach my $otherex (keys %{$exonoverlap->{$exon}}){
		    my $fr1=$exonoverlap->{$exon}{$otherex}{"fractionoverlap1"};
		    my $fr2=$exonoverlap->{$exon}{$otherex}{"fractionoverlap2"};
		    
		    if($fr1 >=0.1 || $fr2 >=0.1){
			
			foreach my $gene (keys %{$exongene->{$otherex}}){
			    my $prefix=substr $gene, 0, 3;
			    
			    if($prefix eq "ENS" && (!exists $monogenes->{$gene})){
				$ensgenes{$gene}=1;
				$ensgenesex{$gene}=1;
			    }
			}
		    }
		}
	    }
	}
	
	my $nbens=keys %ensgenes;
	
	if($nbens>=2){
	    $chimerictx->{$tx}=join(",", keys %ensgenes);
	}
    }
}

##############################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script compares Ensembl and assembled annotations to find readthrough transcripts.\n";
    print "\n";
    
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

##########################################################################
##########################################################################

## parameters 

my %parameters;
$parameters{"pathEnsemblGTF"}="NA";
$parameters{"pathAssembledGTF"}="NA";
$parameters{"pathGeneInfo"}="NA";
$parameters{"monoexonicBiotypes"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathEnsemblGTF", "pathAssembledGTF", "pathGeneInfo", "monoexonicBiotypes", "pathOutput");
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

my %exoncoords;
my %exontx; ## exon - transcript correspondence
my %exongene; ## exon - gene correspondence
my %genetx; ## gene - transcript correspondence
my %txexons; ## transcript-exons correspondence

print "Reading Ensembl GTF...\n";

readGTF($parameters{"pathEnsemblGTF"}, \%exoncoords, \%exontx, \%exongene, \%genetx, \%txexons);
  
my $nbex=keys %exontx;
my $nbgenes=keys %genetx;

print "Found ".$nbgenes." genes and ".$nbex." exons.\n";

print "Done.\n";

##############################################################

print "Reading assembled GTF...\n";

readGTF($parameters{"pathAssembledGTF"}, \%exoncoords, \%exontx, \%exongene, \%genetx, \%txexons);
  
my $nbex2=keys %exontx;
my $nbgenes2=keys %genetx;

print "There are ".$nbgenes2." genes and ".$nbex2." exons after adding assembled annotations.\n";

print "Done.\n";


##############################################################


print "Reading gene info...\n";
my %geneinfo;

readGeneInfo($parameters{"pathGeneInfo"}, \%geneinfo);

print "Done.\n";

print "Extracting genes and exons with monoexonic biotypes...\n";

my @m=split(",", $parameters{"monoexonicBiotypes"});

print join(", ", @m)." are the monoexonic biotypes.\n";

my %monobio;

foreach my $mm (@m){
    $monobio{$mm}=1;
}

my %monoexons;
my %monogenes;

foreach my $gene (keys %genetx){
    if(exists $geneinfo{$gene}){
	
	my $bio=$geneinfo{$gene}{"biotype"};
	
	if(exists $monobio{$bio}){
	    my %thisexon;
	    
	    foreach my $tx (keys %{$genetx{$gene}}){
		foreach my $exid (keys %{$txexons{$tx}}){
		    $thisexon{$exid}=1;
		}
	    }
	    
	    my $nbex=keys %thisexon;

	    if($nbex==1){
		$monogenes{$gene}=1;
		foreach my $exid (keys %thisexon){
		    $monoexons{$exid}=1;
		}
	    }
	    else{
		print "Unusual: ".$gene." is ".$bio." but has ".$nbex." exons.\n";
	    }
	}
    }
}

my $nbmonogenes=keys %monogenes;
my $nbmonoex=keys %monoexons;

print "Found ".$nbmonogenes." and ".$nbmonoex." monoexonic genes and exons.\n";

print "Done.\n";

##############################################################

print "Ordering exons...\n";

my %orderedexons;

orderExons(\%exoncoords, \%orderedexons);

print "Done.\n";

print "Computing exon overlap...\n";
my %exonoverlap;

extractOverlap(\%orderedexons, \%orderedexons, 0, "sense", \%exonoverlap);
 
print "Done.\n";

##############################################################

print "Extracting read-through (chimeric) assembled transcripts...\n";

my %chimerictx;

extractChimericTranscripts(\%txexons, \%exongene, \%monogenes,  \%exonoverlap, \%chimerictx);

my $nbchimerictx=keys %chimerictx;

print "Found ".$nbchimerictx." chimeric assembled transcripts.\n";

print "Done.\n";

##############################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

print $output "TranscriptID\tEnsemblGenes\n";

foreach my $txid (keys %chimerictx){
    my $prefix=substr $txid,0,3;

    ## only newly assembled transcripts

    if($prefix ne "ENS"){ 
	print $output $txid."\t".$chimerictx{$txid}."\n";
    }
}

close($output);

print "Done.\n";

##############################################################
