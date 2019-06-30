use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################
##############################################################

sub readGTF{
    my $pathin=$_[0];
    my $transcripts=$_[1];
    my $exons=$_[2];
    
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

	    my $info=$s[8];
	    my @infoarray=split(";", $info);
	   
	    my $txid=findInfo("transcript_id", \@infoarray);
	    my $exonid=$chr.",".$start.",".$end.",".$strand;
	 
	    if(exists $transcripts->{$txid}){
		push(@{$transcripts->{$txid}{"start"}}, $start);
		push(@{$transcripts->{$txid}{"end"}}, $end);
	    } else{
		$transcripts->{$txid}={"start"=>[$start], "end"=>[$end], "chr"=>$chr, "strand"=>$strand};
	    }

	    if(exists $exons->{$exonid}){
		$exons->{$exonid}{$txid}=1;
	    } else{
		$exons->{$exonid}={$txid=>1};
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
    
    if(@grepres==1){
	my @t=split("\"",$grepres[0]);
	$res=$t[1];
    }
    
    return $res;
}

##############################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script compares de novo transcripts with annotated transcripts.\n";
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
$parameters{"pathKnownGTF"}="NA";
$parameters{"pathAssembledGTF"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathKnownGTF", "pathAssembledGTF", "pathOutput");

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

print "Reading annotations...\n";

my %knowntranscripts;
my %knownexons;

readGTF($parameters{"pathKnownGTF"}, \%knowntranscripts, \%knownexons);

my $nbtx=keys %knowntranscripts;
my $nbex=keys %knownexons;
print "Found ".$nbtx." transcripts and ".$nbex." exons in known GTF.\n";


my %newtranscripts;
my %newexons;

readGTF($parameters{"pathAssembledGTF"}, \%newtranscripts, \%newexons);

my $nbtx=keys %newtranscripts;
my $nbex=keys %newexons;
print "Found ".$nbtx." transcripts and ".$nbex." exons in assembled GTF.\n";

print "Done.\n";

##############################################################

print "Comparing transcripts...\n";

my %commonexons;

foreach my $exonid (keys %knownexons){
    if(exists $newexons{$exonid}){
	foreach my $knowntx (keys %{$knownexons{$exonid}}){
	    foreach my $newtx (keys %{$newexons{$exonid}}){
		if(exists $commonexons{$knowntx}){
		    if(exists $commonexons{$knowntx}{$newtx}){
			$commonexons{$knowntx}{$newtx}{$exonid}=1;
		    }
		    else{
			$commonexons{$knowntx}{$newtx}={$exonid=>1};
		    }
		}
		else{
		    $commonexons{$knowntx}={$newtx=>{$exonid=>1}};
		}
	    }
	}
    }
}
print "Done.\n";

##############################################################

print "Writing output for synonyms...\n";

open(my $output, ">".$parameters{"pathOutput"});

print $output "DeNovoID\tEnsemblID\n";

my $nbsyn=0;

foreach my $knowntx (keys %commonexons){
    my $nbknown=@{$knowntranscripts{$knowntx}{"start"}};
    
    foreach my $newtx (keys %{$commonexons{$knowntx}}){
	my $nbnew=@{$newtranscripts{$newtx}{"start"}};

	my $nbcommon=keys %{$commonexons{$knowntx}{$newtx}};

	if($nbnew==$nbknown && $nbcommon==$nbnew){
	    print $output $newtx."\t".$knowntx."\n";
	}
    }
}

close($output);

print "Done.\n";

##############################################################
