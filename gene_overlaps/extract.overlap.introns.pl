use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

###################################################################
###################################################################

sub readGTF{
    my $pathin=$_[0];
    my $forbidden=$_[1];
    my $genecoords=$_[2];
    
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

		if(exists $forbidden->{$tx}){
		    print "Discarding transcript ".$tx."\n";
		}
		else{
		    if($gene eq "NA"){
			print "Weird! cannot find gene id for ".$line."\n";
			exit(1);
		    }
		    
		    if(exists $genecoords->{$gene}){
			if($start < $genecoords->{$gene}{"start"}){
			    $genecoords->{$gene}{"start"}=$start;
			}
			
			if($end > $genecoords->{$gene}{"end"}){
			    $genecoords->{$gene}{"end"}=$end;
			}
			
		    } else{
			$genecoords->{$gene}={"start"=>$start, "end"=>$end, "chr"=>$chr, "strand"=>$strand};
		    }
		}
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

##############################################################

sub readIntronCoords{
    my $pathin=$_[0];
    my $introns=$_[1];

    open(my $input, $pathin);
    my %header;
    my $line=<$input>;
    chomp $line;
    
    my @s=split("\t", $line);
    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }

    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $chr=$s[$header{"Chr"}];
	my $start=$s[$header{"Start"}];
	my $end=$s[$header{"End"}];
	my $strand=$s[$header{"Strand"}];

	my $intid=$chr.",".$start.",".$end.",".$strand;

	$introns->{$intid}={"chr"=>$chr, "start"=>$start, "end"=>$end, "strand"=>$strand};
	
	$line=<$input>;
    }
    
    close($input);
}

##############################################################

sub orderCoords{
    my $coords=$_[0];
    my $ordered=$_[1];
    
    my %hashpos;

    foreach my $exid (keys %{$coords}){
	my $chr=$coords->{$exid}{"chr"};
	my $strand=$coords->{$exid}{"strand"};
	my $start=$coords->{$exid}{"start"};
	my $end=$coords->{$exid}{"end"};

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

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script computes overlaps between two sets of genes. \n";
    print "\n";
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

##################################################################
##################################################################
