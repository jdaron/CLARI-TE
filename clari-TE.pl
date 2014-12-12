#!/usr/local/bin/perl

#~ Licence
#~ =======
#~ Copyright or © or Copr. Josquin Daron INRA-GDEC 01/09/2014
#~ 
#~ email: josquin.daron@clermont.inra.fr
#~ 
#~ This software is a computer program whose purpose is to predict Transposable 
#~ Elements (TEs) in complexe genome such as wheat. The program correct raw similarity
#~ search result from RepeatMasker.
#~ 
#~ This software is governed by the CeCILL license under French law and
#~ abiding by the rules of distribution of free software.  You can  use, 
#~ modify and/ or redistribute the software under the terms of the CeCILL
#~ license as circulated by CEA, CNRS and INRIA at the following URL
#~ "http://www.cecill.info". 
#~ 
#~ As a counterpart to the access to the source code and  rights to copy,
#~ modify and redistribute granted by the license, users are provided only
#~ with a limited warranty  and the software's author,  the holder of the
#~ economic rights,  and the successive licensors  have only  limited
#~ liability. 
#~ 
#~ In this respect, the user's attention is drawn to the risks associated
#~ with loading,  using,  modifying and/or developing or reproducing the
#~ software by the user in light of its specific status of free software,
#~ that may mean  that it is complicated to manipulate,  and  that  also
#~ therefore means  that it is reserved for developers  and  experienced
#~ professionals having in-depth computer knowledge. Users are therefore
#~ encouraged to load and test the software's suitability as regards their
#~ requirements in conditions enabling the security of their systems and/or 
#~ data to be ensured and,  more generally, to use and operate it in the 
#~ same conditions as regards security. 
#~ 
#~ The fact that you are presently reading this means that you have had
#~ knowledge of the CeCILL license and that you accept its terms.

use strict ;
use warnings ;
use Bio::SeqIO ;
use diagnostics ;
use Bio::SeqFeature::Generic ;
use Bio::Tools::GFF ;
use Data::Dumper ;
use File::Basename ;
use Statistics::Descriptive ;
use Getopt::Long ;
use Bio::Location::Simple ;

my $VERSION = "1" ;
my $lastmodif = "01/02/2014" ;
my $help ;
my $annotLTR_File ;
my $classiFile ;
my $optFasta ;
my $outputFormat = "embl" ;
my $geneEmblFile  ;
my $geneFormat = "triAnnot" ;
my $dir = "" ;
my $succefully = "--> failure is not an option\n" ;
my $verbosity = 3 ;

=head
CLARI-TE Option: 

-h|--help:             print this help
-LTR:                  tabulation file of the LTR annotations
-classi:               corresponding file sequence id to family id
-fasta:                fasta file of the annotated sequence
-gene:                 embl file of the gene annotation
-dir:                  directory name
-v:                    verbosity (3,4)
=cut

&GetOptions ( "h|help"      => \$help,
			  "LTR:s"       => \$annotLTR_File,
			  "classi:s"    => \$classiFile,
			  "fasta:s"     => \$optFasta,
			  "gene:s"      => \$geneEmblFile,
			  "dir:s"       => \$dir,
			  "v:i"         => \$verbosity);

$help    and &help ;
@ARGV     or print STDERR "*** No argument file (xm only) ***\n" and &help;
my $source_tag = "post_RM" ;
my $kp = 0 ;
my $tmpK = 0 ;

&main(\@ARGV);
exit(0) ;

sub main {
	my $self = {};
	bless $self;
	$self->{xmInputFiles} = shift;
	$self->getClassification() ;
	$self->getFasta() ;
	$self->getLTRposition() ;
	foreach (@{$self->{xmInputFiles}}){
		$self->{_inputFile} = $_;
		$self->getTeAnnotation() ;
		$self->getGeneAnnotation() ;
		$self->sortFeat() ;
		$self->step1ResolveOverlap() ;
		$self->step2MergeFeatures() ;
		$self->step3JoinFeatures() ;
	}
	exit(0);
}

sub getClassification {
	my $self = shift;
	if ( not defined ($classiFile)) { die "*** Could not find file containing the classification of the library (option -classi)  ***\n" ; }
	open (LIST, $classiFile) or die "*** Could not open file containing the classification of the library : \"$classiFile\" ***\n" ;
	while(<LIST>) {
		chomp $_ ;
		if ($_=~/^#/) {next ;}
		my @tab = split(" ", $_) ;
		unless (scalar(@tab) == 2)  { die ("*** Problems of format in \"$classiFile\" : wrong number of columns, need 2 columns ***\n") ; } 
		$self->{_classification}->{$tab[0]}->{fam} = $tab[1] ;
	}
	close LIST ;
}

sub getFasta { # Read FASTA file
	my $self = shift;
	$self->{fasta}->{name} = $optFasta ;
	my $seqio_obj = Bio::SeqIO->new(-file => $self->{fasta}->{name}, -format => "FASTA");
	$self->{fasta}->{seq} = $seqio_obj->next_seq->seq ;
}

sub getLTRposition{ # Read tabe file of LTR annotation
	my $self = shift;
	if (not defined ($annotLTR_File) ) {die ("*** Could not find file of position of LTR (option -LTR)***\n") ; }
	open(LTR, $annotLTR_File) or die ("*** Can't open file of position of LTR \"$annotLTR_File\" ***\n")  ;
	while(<LTR>) {
		chomp $_ ;
		if ($_=~/^#/) { next; }
		
		my @tab = split("\t",$_) ;
		unless (scalar(@tab) == 5) { die ("*** Problems of format in \"$annotLTR_File\" : wrong number of columns, need 5 columns ***\n"); }
		unless ($tab[2] eq "LTR5" or $tab[2] eq "LTR3") { die ("*** Problems of format in \"$annotLTR_File\" : third column \"$tab[2]\" should be \"LTR 5'\" or \"LTR 3'\" ***\n"); }
		$self->{_classification}->{$tab[0]}->{$tab[2]}->{start} = $tab[3] ; 
		$self->{_classification}->{$tab[0]}->{$tab[2]}->{end} = $tab[4] ;
	}
}

sub getTeAnnotation{ # Read TE annotation, only xm format, but other type of format (TEannot,...) could be added !
	my $self = shift;
	unless (-e $self->{_inputFile}) { print STDERR ("*** Could not find file in argument (xm) : \"$self->{_inputFile}\" ***\n") and die ;  }
	my @annotTE ;
	if ($self->{_inputFile}=~/^(.*).out.xm$/) {
		$self->readXM();
	} else { print STDERR ("*** can not recognize file (or file's extention) : \"$self->{_inputFile}\" ; just xm are accepted\n" ) and die; }
	my $outfile = $1 ;
	$self->{output}->{prefixName} = basename($outfile) ;
	foreach my $obj (@{$self->{inputFile}->{xm}}) { 	# group feature in a tmp_hash
		if ($obj->source_tag eq "xm") { my $annotType = "RM" ; $self->{tmp_feat}->{++$tmpK} = $obj ; } # part dedicate to RepeatMasker
	}
}

sub getGeneAnnotation{ # put Gene annotation in a format legible to clariTE.pl, only triAnnot format is accepted, but other type of format could be added !
	my $self = shift;
	my $totalEmblObj = 0 ;
	if (defined $geneEmblFile) {
		$self->readEMBL($geneEmblFile) ;
		$totalEmblObj = scalar(@{$self->{inputFile}->{gene}}) ;
		foreach my $obj (@{$self->{inputFile}->{gene}}) {
			my $info ;
			my $ID ;
			if ($geneFormat eq "triAnnot") {
				unless ($obj->primary_tag eq "CDS" ) { next ;  }
				$info->{targetFamID} = "exon," ;
				$info->{targetVarID} = "exon," ;
				$info->{copieID} = "gene" ;
				$info->{targetLgth} = $obj->length ;
				$info->{targetStart} = 1 ;
				$info->{targetEnd} = $obj->length ;
				$obj->{info} =  $info ;
				if ($obj->has_tag("locus_tag")) {$ID = ($obj ->each_tag_value("locus_tag"))[0] ; }
				my $note = "exon ".$obj->length."bp 1..".$obj->length ;
				$obj ->add_tag_value ( "id", $ID ) ;
				$obj ->add_tag_value ( "post", $note ) ;
				$self->{tmp_feat}->{++$tmpK} = $obj ;
			}
		}
	} else { 
		print STDERR ("*** WARNING :No gene annotation file give as option -gene ***\n") ; 
	}
}

sub sortFeat { # sort all features from the TE annotation and Gene annotation
	my $self = shift;
	my $k = 0 ;
	foreach my $matchID ( sort { $self->{tmp_feat}->{$a}->start <=> $self->{tmp_feat}->{$b}->start || $self->{tmp_feat}->{$a}->end <=> $self->{tmp_feat}->{$b}->end } keys %{$self->{tmp_feat}}) {
		$k++ ;
		$self->{feat}->{$k} = $self->{tmp_feat}->{$matchID} ;
	}
	if (scalar (keys %{$self->{feat}}) < 2 and scalar (keys %{$self->{feat}}) > 0) {  # not enough feature to continue Clari-TE
		$self->{_print} = $self->{feat} ;
		$self->{_tag} =  "_annoTE." ;
		$self->finishClariTE() ;
	}
	if ($verbosity == 4 ) {
		$self->{_print} = $self->{feat} ;
		$self->{_tag} =  "_outRM." ;
		$self->finishClariTE() ;
	}
}

sub finishClariTE{ # print features by creating a EMBL file
	my $self = shift ;
	delete $self->{final} ;
	foreach my $matchID (sort { $self->{_print}->{$a}->start <=> $self->{_print}->{$b}->start } keys %{$self->{_print}} ) {
		push (@{$self->{final}}, $self->{_print}->{$matchID}) ;
	}
	my $outputName = $self->{output}->{prefixName}.$self->{_tag}.$outputFormat ;
	if ($outputFormat eq "embl") {
		&createEmbl (\@{$self->{final}}, $outputName, $self->{fasta}->{seq} ) ;
	}
	if ($self->{_tag} eq "_annoTE.") {
		print $succefully ;
		exit(0) ;
	}
}

sub step1ResolveOverlap { # STEP 1
	print STDERR ("--> step overlaping feature\n") ;
	my $self = shift ;
	$self->overlaping() ;

	if (scalar (keys %{$self->{noOverlapingFeat}}) < 3) { # not enough feature to continue Clari-TE
		$self->{_print} = $self->{noOverlapingFeat} ;
		$self->{_tag} =  "_annoTE." ;
		$self->finishClariTE() ;
	}

	if ($verbosity == 4 ) {
		print STDERR ("--> check if still overlaping feature\n") ;
		for (my $i = 2 ; $i < scalar ( keys (%{$self->{noOverlapingFeat}} )) ; $i ++ ) {
			if ( $self->{noOverlapingFeat}->{$i-1}->end > $self->{noOverlapingFeat}->{$i}->start ) { # check if feature still overlapping 
				print STDERR ("WARNING : Overlaping Features : \n".$self->{noOverlapingFeat}->{$i-1}->start."\t".$self->{noOverlapingFeat}->{$i-1}->end."\t".$self->{noOverlapingFeat}->{$i}->start."\t".$self->{noOverlapingFeat}->{$i}->end."\n") ;
			}
		}
		$self->{_print} = $self->{noOverlapingFeat} ;
		$self->{_tag} =  "_resolveOverlap." ;
		$self->finishClariTE() ;
	}

	# step filtering small features
	my $tmpK = 0 ;
	for (my $i = 1 ; $i <= scalar ( keys (%{$self->{noOverlapingFeat}} )) ; $i ++ ) { 
		if ($i == scalar (keys (%{$self->{noOverlapingFeat}}))) { $tmpK++ ; $self->{noOverlapingFeat_noSmallFeat}->{$tmpK} = $self->{noOverlapingFeat}->{$i} ; last ; } # last round of the loop
		my $cutoff = 50 ;
		$self->{noOverlapingFeat}->{$i}->{compo}->{$self->{noOverlapingFeat}->{$i}->{info}->{targetVarID}}->{lgth} = $self->{noOverlapingFeat}->{$i}->length ;
		$self->{noOverlapingFeat}->{$i}->{compo}->{$self->{noOverlapingFeat}->{$i}->{info}->{targetVarID}}->{conslgth} = $self->{noOverlapingFeat}->{$i}->{info}->{targetLgth} ;
		my $test = &filterOutSmallFeat ($i , $self->{noOverlapingFeat}, $cutoff) ;
		if ($test == 0) {
			$tmpK++ ;
			$self->{noOverlapingFeat_noSmallFeat}->{$tmpK} = $self->{noOverlapingFeat}->{$i} ;
		}
	}
	if (scalar (keys %{$self->{noOverlapingFeat_noSmallFeat}}) < 3) { # not enough feature to continue Clari-TE
		$self->{_print} = $self->{noOverlapingFeat_noSmallFeat} ;
		$self->{_tag} =  "_annoTE." ;
		$self->finishClariTE() ;
	} 
}

sub step2MergeFeatures{ # STEP 2
	print STDERR ("--> step merge collinear feature \n") ;
	my $self = shift ;
	if ($verbosity == 4) { print STDERR ("--> resolveLTR\n\tTEid\tstart\tstrand\tend\ttargetStart\ttargetEnd\n") ;}
	my $kFeat = 0 ;
	my $flag = 0 ;

	for (my $i = 1 ; $i <= scalar (keys (%{$self->{noOverlapingFeat_noSmallFeat}})) ; $i++ ) {
		$kFeat++ ;
		if ($i == 1 ) { $self->{mergeFeat}->{$kFeat} = $self->{noOverlapingFeat_noSmallFeat}->{$i} ; next ; } # first round of the loop
		my ($prev, $cur, $next) ; 
		my $newFeat = 0 ;
		
		TOP:
		if ( $flag == 0 ) { 
			$prev = $self->{mergeFeat}->{$kFeat-1} ;
			$cur = $self->{noOverlapingFeat_noSmallFeat}->{$i} ; 
			$next = $self->{noOverlapingFeat_noSmallFeat}->{$i+1} ; 
		}
		elsif ($flag == 1 ) { 
			$prev = $self->{mergeFeat}->{$kFeat-1} ;
			$cur = $newFeat ; 
			$next = $self->{noOverlapingFeat_noSmallFeat}->{$i+1} ; 
		}
		
		if ($i == scalar (keys (%{$self->{noOverlapingFeat_noSmallFeat}}))) { # last round of the loop
			my $feat_Col = &mergeColFeat ($prev, $cur, 0, $self) ;
			if (defined ($feat_Col->{cur})) { $self->{mergeFeat}->{$kFeat} = $feat_Col->{cur} ; }
			else { $kFeat = $kFeat - 1 ; $self->{mergeFeat}->{$kFeat} = $feat_Col->{newFeat} ; }
			last ; 
		}
		$flag = 0 ;
		
		my $feat_Col = &mergeColFeat ($prev, $cur, $next, $self) ;
		if (defined ($feat_Col->{cur})) { # not merge of features
			my $wrongFeat = &filterOutWongMatch ($prev, $cur, $next) ;
			if ($wrongFeat == 0 ) { # $cur is conserved 
				$self->{mergeFeat}->{$kFeat} = $feat_Col->{cur} ;
			}
			else { # cur match a removed, merge feature prev with next
				$i = $i + 1 ;
				$kFeat = $kFeat-1 ;
				$newFeat = $wrongFeat ;
				if ($kFeat == 1 ) { $self->{mergeFeat}->{$kFeat} = $newFeat ; next ; }
				$flag = 1 ;
				goto (TOP) ;
			}
		}
		else { # merge features prev and cur
			$kFeat = $kFeat - 1 ;
			$newFeat = $feat_Col->{newFeat} ;
			if ($kFeat == 1 ) { $self->{mergeFeat}->{$kFeat} = $newFeat ; next ; }
			$flag = 1 ;
			goto (TOP) ;
		}
	}
	
	if (scalar (keys %{$self->{mergeFeat}}) < 3) {
		$self->{_print} = $self->{mergeFeat} ;
		$self->{_tag} =  "_annoTE." ;
		$self->finishClariTE() ;
	} 

	if ($verbosity == 4 ) {
		$self->{_print} = $self->{mergeFeat} ;
		$self->{_tag} =  "_mergeFeat." ;
		$self->finishClariTE() ;
	}
}

sub step3JoinFeatures { # STEP 3
	print STDERR ("--> step Join \n") ;
	my $self = shift ;
	if ($verbosity == 4) { print STDERR (scalar (keys (%{$self->{mergeFeat}}))." subfeature after merge step\n") ; }
	my $round = 0 ;
	my $j = 1 ;
	JOIN:
	my $joinSize ;
	if ($round == 0) { $joinSize = 0 ; }
	else { $joinSize = scalar(@{$self->{joinFeat}}) ;}
	my $nbreJoin = 0 ;
	for (my $i = 2 ; $i + $j < scalar (keys (%{$self->{mergeFeat}})) ; $i++ ) {
		if ($self->{mergeFeat}->{$i}->{info}->{targetFamID} eq "exon" and $self->{mergeFeat}->{$i-1}->{info}->{targetFamID} eq "exon") {
			if ($self->{mergeFeat}->{$i}->has_tag("locus_tag") and $self->{mergeFeat}->{$i-1}->has_tag('locus_tag')) {
				if (($self->{mergeFeat}->{$i}->each_tag_value('locus_tag'))[0] eq ($self->{mergeFeat}->{$i-1}->each_tag_value('locus_tag'))[0] ) {
					my $feature = &createParentFeature ( \$self->{mergeFeat}->{$i-1}, \$self->{mergeFeat}->{$i}) ;
					$feature->add_tag_value ("locus_tag", ($self->{mergeFeat}->{$i}->each_tag_value('locus_tag'))[0] ) ;
					push (@{$self->{joinFeat}}, $self->{mergeFeat}->{$i-1}) ;
					push (@{$self->{joinFeat}}, $self->{mergeFeat}->{$i}) ;
					$self->{mergeFeat}->{$i} = $feature ;
					delete ( $self->{mergeFeat}->{$i-1} ) ;
					$nbreJoin++ ;
					next ;
				}
			}
		}
		my $test =0 ;
		if ($self->{mergeFeat}->{$i}->{info}->{targetFamID} eq "gap") {
			$test = &joinFeat($i, $j, $self->{mergeFeat}, 2000 ) ;
		}
		elsif ($self->{mergeFeat}->{$i}->{info}->{targetStart} < 50 and $self->{mergeFeat}->{$i}->{info}->{targetEnd} + 50 > $self->{mergeFeat}->{$i}->{info}->{targetLgth}) { # feature $i is complete (join feature trough complete TE)
			$test = &joinFeat($i, $j, $self->{mergeFeat}, 2000 ) ;
		}
		elsif ($self->{mergeFeat}->{$i}->primary_tag eq "match_part") {
			$test = &joinFeat($i, $j, $self->{mergeFeat}, 1000 ) ;
		}
		else { $test = &joinFeat($i, $j, $self->{mergeFeat}, 500 ) ; }
		if ($test == 1) {
			my ($feature) = &createParentFeature ( \$self->{mergeFeat}->{$i-1}, \$self->{mergeFeat}->{$i+$j}) ;
			for (my $k = $i-1 ; $k <= $i+$j ; $k++) {
				unless ( ($self->{mergeFeat}->{$k}->each_tag_value('id'))[0] eq ($feature->each_tag_value('id'))[0]) { push (@{$self->{joinFeat}}, $self->{mergeFeat}->{$k}) ; }
				delete ($self->{mergeFeat}->{$k}) ;
			}
			$self->{mergeFeat}->{$i+$j} = $feature ;
			$i = $i + $j ;
			$nbreJoin += $test ;
		}
	}
	$round++ ;
	my $k = 0 ;
	my $size = scalar (keys (%{$self->{mergeFeat}})) ;
	foreach my $key ( sort { $self->{mergeFeat}->{$a}->start <=> $self->{mergeFeat}->{$b}->start } keys %{$self->{mergeFeat}}) {
		$k++ ;
		$self->{mergeFeat}->{$k} = $self->{mergeFeat}->{$key} ;
	}
	foreach my $key ( sort { $self->{mergeFeat}->{$a}->start <=> $self->{mergeFeat}->{$b}->start } keys %{$self->{mergeFeat}}) {
		unless ($key <= $size) { delete $self->{mergeFeat}->{$key} ; }
	}
	if ($joinSize != scalar (@{$self->{joinFeat}})) { goto(JOIN) ; } 
	if ($j < 10 ) { $j++ ; goto(JOIN) ; } 
	
	$k = 0 ;
	$size = scalar (keys (%{$self->{mergeFeat}})) ;
	foreach my $obj (@{$self->{joinFeat}}) {
		$k++ ;
		$self->{mergeFeat}->{$k+$size} = $obj ;
	}
	
	$self->{_print} = $self->{mergeFeat} ;
	$self->{_tag} =  "_annoTE." ;
	$self->finishClariTE() ;
}

#===== SUB Read and Create files =====

sub readXM {
	my $self = shift ;
	open (XM, $self->{_inputFile}) or die print STDERR ("*** Can't open file \"$self->{_inputFile}\" ***\n") ;
	my $id = 0 ;
	while (<XM>){
		my @col = split (/\s+/, $_) ;
		$id++ ;
		my $prefixFileName = $col[4] ;
		my $query_name     = $col[4] ;
		my $rep_score      = $col[0] ;
		my $q_match_start  = $col[5] ;
		my $q_match_end    = $col[6] ;
		my $q_match_strand = $col[8] ;
		my $rep_name       = $col[9] ;
		$rep_name =~/^(.*)#Unknown$/i and $rep_name = $1 ;
		my $h_match_start ;
		my $h_match_end ;
		my $h_unmasked ;
		my $ft_color ;
		my $rep_strand ;
		my $info ;
		if ($q_match_strand eq "+") {
			$rep_strand = 1 ;
			$h_match_start  = $col[10] ;
			$h_match_end    = $col[11] ;
			$h_unmasked     = $col[12] ;
		}
		else {
			$rep_strand = -1 ;
			$h_match_start  = $col[12] ;
			$h_match_end    = $col[11] ;
			$h_unmasked     = $col[10] ;
		}
		$h_unmasked =~s/\(//g ;
		$h_unmasked =~s/\)//g ;
		my $hitlen = $h_match_end+$h_unmasked ;
		
		if (defined($self->{_classification}->{$rep_name}->{fam})){
			if ($self->{_classification}->{$rep_name}->{fam}=~/([A-Z]*_fam[n|c][0-9]*)(\.[0-9]*)/) { 
				$info->{targetFamID} = $1 ;
				$info->{targetVarID} = $1.$2 ;
			}
			elsif ($self->{_classification}->{$rep_name}->{fam}=~/([A-Z]*_fam[n|c][0-9]*)/) { 
				$info->{targetFamID} = $1 ;
				$info->{targetVarID} = $1 ;
			}
			else { print STDERR ("can not find \"$self->{_classification}->{$rep_name}->{fam}\" in classification !\n") ; next ; }

			$info->{copieID} = $rep_name ;
			$info->{targetLgth} = $hitlen ;
			$info->{targetStart} = $h_match_start ;
			$info->{targetEnd} = $h_match_end ;
			
			my $feat = new Bio::SeqFeature::Generic ( -primary      => 'repeat_region',
													  -start        => $q_match_start,
													  -end          => $q_match_end,
													  -strand       => $rep_strand,
													  -source_tag   => 'xm') ;
			my $post ;
			if (defined ($self->{_classification}->{$rep_name})){ $post = $self->{_classification}->{$rep_name}->{fam}." ".$hitlen."bp ".$h_match_start."..".$h_match_end ; }
			else { $post = $rep_name." ".$hitlen."bp ".$h_match_start."..".$h_match_end ; }
			$feat ->add_tag_value ("post", $post) ;
			$feat->{info} =  $info ;
			$feat ->add_tag_value ("id", $id."_".$prefixFileName) ;
			push (@{$self->{inputFile}->{xm}}, $feat) ;
		}
		else { print STDERR ("can not find \"$rep_name\" in classification !\n") ; next ; }
	}
}

sub readEMBL {
	my $self = shift ;
	my $file = $_[0] ;
	-e $file or die "*** file \"$file\" not found ***\n";
	my $seqIO_obj = Bio::SeqIO->new( -format => "EMBL", -file => $file) ;
	while ( my $seq_obj = $seqIO_obj->next_seq ) {
		my @features = $seq_obj->all_SeqFeatures() ;
		foreach my $feature (@features) {
			my $pritag = $feature->primary_tag ;
			my $strand = $feature->strand ;
			my $location = $feature->location ;
			my $seqID = $feature->seq_id ;

			
			my @locs = sort {$a->start <=> $b->start} $location->each_Location ; # listes des localisations (start end) des subfeatures pour chaque features 
			foreach my $loc (@locs) {
				my $start = $loc->start ;
				my $end = $loc->end ;
				my $newFeat = Bio::SeqFeature::Generic -> new ( -seq_id      => $seqID,
																-source_tag  => $source_tag,
																-primary_tag => $pritag,
																-start       => $start,
																-end         => $end,
																-strand      => $strand
																 ); 
				if ($feature->has_tag ("locus_tag")) {
					my $locus_tag = ($feature ->each_tag_value("locus_tag"))[0] ;
					$newFeat ->add_tag_value ( "locus_tag", $locus_tag ) ;
				}
				elsif ($feature->has_tag ("ID")) {
					my $locus_tag = ($feature ->each_tag_value("ID"))[0] ;
					$newFeat ->add_tag_value ( "locus_tag", $locus_tag ) ;
				}
				$newFeat ->add_tag_value ( "exonNbre", scalar (@locs) ) ;
				push (@{$self->{inputFile}->{gene}}, $newFeat) ;
			}
		}
	}
}
 
sub createEmbl {
	my ($final, $outfile, $seq) = @_ ;
	my $k = 0 ;
	my $featPart ;
	my @feat ;
	my $seq_id ;
	my $match_part ;
	foreach $_ (@{$final}) {
		$seq_id = $_->seq_id ;
		my $pritag = $_->primary_tag ;

		if ($pritag eq "match_part") {
			my $parent = ($_->each_tag_value('id'))[0] ;
			$match_part->{$parent}->{post} = ($_->each_tag_value('post'))[0] ;
			if (defined ($_->{info})) { $match_part->{$parent}->{info} = $_->{info} ; }
			if (defined ($_->{range})) {
				my $range ;
				my $kr =0 ;
				foreach my $start (sort { $a <=> $b } keys %{$_->{range}}) {
					$kr++ ;
					$range .= $_->{range}->{$start} ;
					unless ($kr == scalar(keys(%{$_->{range}}))) { $range .= "," ; }
				}
				$match_part->{$parent}->{range} = $range ;
			}
			next ;
		}
		if ( $_->has_tag("parent") ) {
			my $parent = ($_->each_tag_value('parent'))[0] ;
			$featPart->{$parent}->{$pritag}->{++$k} = $_ ;
		}
		else {
			if ($_->primary_tag eq "repeat_region") {
				my $compo = "" ;
				my $totalPercent = 0 ;
				foreach my $var (keys %{$_->{compo}}) {
					my $percent = sprintf "%.2f",($_->{compo}->{$var}->{lgth}/$_->length)*100 ; 
					$compo .= $var." ". $percent." " ;
					$totalPercent += $percent ;
				}
				if (100 - $totalPercent != 0) { $compo .= "no_match ". sprintf "%.2f",(100 - $totalPercent)." " ; }
				if ($_->has_tag('compo')) { $_->remove_tag('compo') ; }
				$_->add_tag_value('compo', "$compo" );
				if (defined ($_->{soloLTR})) {
					if ($_->has_tag('soloLTR')) { $_->remove_tag('soloLTR') ; }
					$_->add_tag_value('soloLTR', 'soloLTR');
				}
			}
			push (@feat, $_ ) ; 
		}
	}
	
	foreach my $parent (sort keys %{$featPart}) {
		PRITAG:
		foreach my $pritag (sort keys %{$featPart->{$parent}}) {
			my $splitlocation = Bio::Location::Split->new ;
			my $source ;
			my $strand ;
			my $k2 = 0 ;
			my @coord ;
			my $ID ;
			my $note ;
			my $parentTag ;
			my $label ;
			my $range ;
			if ($pritag eq "repeat_region") { $range = $match_part->{$parent}->{range} ; }
			my $compo ;
			my $status ;
			my $lgthFeat = 0 ;
			my $composition ;
			my $info ;
			if ($pritag eq "repeat_region") {$info = $match_part->{$parent}->{info} ; }
			
			foreach my $k ( sort{$a<=>$b} keys %{$featPart->{$parent}->{$pritag}} ) {
				$ID = ($featPart->{$parent}->{$pritag}->{$k}->each_tag_value('id'))[0] ;
				$parentTag = ($featPart->{$parent}->{$pritag}->{$k}->each_tag_value('parent'))[0] ;
				$note = $match_part->{$parent}->{post} ;
				$label = $match_part->{$parent}->{label} ;
				
				if ($pritag eq "repeat_region") {
					foreach my $var (keys %{$featPart->{$parent}->{$pritag}->{$k}->{compo}}) {
						if (not defined ($composition->{$var})) { $composition->{$var}->{lgth} = 0 ; $composition->{$var}->{conslgth} = 0 ; }
						$composition->{$var}->{lgth} += $featPart->{$parent}->{$pritag}->{$k}->{compo}->{$var}->{lgth} ;
						$composition->{$var}->{conslgth} = $featPart->{$parent}->{$pritag}->{$k}->{compo}->{$var}->{conslgth} ;
					}
				}
				
				$k2 += 1 ;
				if ($k2 == 1) { $strand = $featPart->{$parent}->{$pritag}->{$k}->strand ; }
				else {
					unless ( $featPart->{$parent}->{$pritag}->{$k}->strand eq $strand ) {
						print STDERR "*** different \"strands\" for JOIN feature with Parent=$parent\n" ;
						next PRITAG ;
					}
				}
				$lgthFeat += $featPart->{$parent}->{$pritag}->{$k}->length ;
				my $start = $featPart->{$parent}->{$pritag}->{$k}->start ;
				my $end   = $featPart->{$parent}->{$pritag}->{$k}->end ;
				my $location_obj = Bio::Location::Simple->new ( '-start' => $start, '-end' => $end, '-strand' => $strand ) ;
				$splitlocation->add_sub_Location( $location_obj ) ;
			} # END OF FOREACH $k
			#-------------------------------------------------------------------------------------
			# This part is used to sort the sub features according to their start position
			my @sublocs = $splitlocation->sub_Location ;
				
			$sublocs[0]->strand == 1  and @sublocs = sort {$a->start <=> $b->start} @sublocs ;
			$sublocs[0]->strand == -1 and @sublocs = sort {$b->start <=> $a->start} @sublocs ;
			my $splitlocation_sort = Bio::Location::Split->new ;
			$splitlocation_sort->add_sub_Location(@sublocs) ;
			#-------------------------------------------------------------------------------------
			my $emblfeat = Bio::SeqFeature::Generic ->new( -seq_id      => $seq_id,
														   -source_tag  => $source,
														   -primary_tag => $pritag,
														   -location    => $splitlocation_sort,
														  ) ;
														  
			if ($pritag eq "repeat_region") {
				my $compo = "" ;
				my $totalPercent = 0 ;
				foreach my $var (keys %{$composition}) {
					my $percent = sprintf "%.2f",($composition->{$var}->{lgth}/$lgthFeat)*100 ; 
					$compo .= $var." ". $percent." " ;
					$totalPercent += $percent ;
				}
				if (100 - $totalPercent != 0) { $compo .= "no_match ". sprintf "%.2f",(100 - $totalPercent)." " ; }
				$emblfeat->add_tag_value('compo', "$compo" );
			}
			$emblfeat->add_tag_value('id', "$ID" );
			$emblfeat->add_tag_value('post', "$note" );
			$emblfeat->add_tag_value('parent', "$parent" );
			$emblfeat->{info} = $info ;
			if (defined ($range)) { $emblfeat->add_tag_value('range', "$range" ); }
			push (@feat, $emblfeat) ;
		}
	}
	@feat = sort {$a->start <=> $b->start} @feat ;

	foreach $_ (@feat) {
		if ($_->primary_tag eq "repeat_region") {
			if ($_->has_tag('copie')) { $_->remove_tag('copie') ; }
			$_->add_tag_value('copie', $_->{info}->{copieID}) ;
			if ($_->has_tag('status')) { $_->remove_tag('status') ; }
			if ($_->has_tag('parent')) {
				my @subLocation = $_->location->each_Location ;
				my $lgthFeat ;
				foreach my $subloc (@subLocation) {
					$lgthFeat += $subloc->length ;
				}
				if ((abs($lgthFeat/($_->{info}->{targetLgth})) > 0.90)) {
					$_->add_tag_value('status', 'complete');
				}
				elsif ($_->{info}->{targetStart} < 50 and $_->{info}->{targetEnd} + 50 > $_->{info}->{targetLgth} and (abs($lgthFeat/$_->{info}->{targetLgth}) > 0.5) ) {
					$_->add_tag_value('status', 'complete');
				}
				else { $_->add_tag_value('status', 'fragmented') ; }
			}
			else {
				if ((abs($_->length/$_->{info}->{targetLgth}) > 0.90)) {
					$_->add_tag_value('status', 'complete');
				}
				elsif ($_->{info}->{targetStart} < 50 and $_->{info}->{targetEnd} + 50 > $_->{info}->{targetLgth} and(abs($_->length/$_->{info}->{targetLgth}) > 0.5)) {
					$_->add_tag_value('status', 'complete');
				}
				else { $_->add_tag_value('status', 'fragmented') ; }
			}
		}
	}
	
	my $seq_obj  = Bio::Seq->new ( '-id' =>  "unknown", '-seq' => $seq ) ;
	$seq_obj->add_SeqFeature(@feat);
	my $seqio_obj = Bio::SeqIO->new ( -file => ">$dir$outfile" , -format => "EMBL", -verbose => -1) ;
	$seqio_obj->write_seq ($seq_obj);
	print STDERR ("--> ". $dir."".$outfile . " embl file created\n");
	return 1 ;
}

#===== SUB Step 1 =====
sub resolveOverlap {
	my ($oldFeat, $currentFeat) = @_ ;
	if (not defined ($oldFeat)) {
		print STDERR ( "*** old obj is undefined value\n" ) and die ;
		if (not defined ($currentFeat)) {
			print STDERR ( "*** current obj is undefined value" ) and die ;
		}
	}
	my $feature ;
	if ( $oldFeat->{info}->{copieID} eq "gene" ) { # mite ou feature complete
		if ($currentFeat->strand eq $oldFeat->strand ) {
			$feature = &newSt_sameStrand ($currentFeat, $oldFeat) ;
		}
		elsif ($currentFeat->strand eq 1 and $oldFeat->strand eq -1) {
			$feature = &newSt_keepReverse ($currentFeat, $oldFeat) ;
		}
		elsif ($currentFeat->strand eq -1 and $oldFeat->strand eq 1) {
			$feature = &newSt_keepReverse ($currentFeat, $oldFeat) ;
		}
		return ($oldFeat, $feature) ;
	}
	elsif ( $currentFeat->{info}->{copieID} eq "gene" ) { # mite ou feature complete
		if ($currentFeat->strand eq $oldFeat->strand ) {
			$feature = &newSp_sameStrand ($oldFeat, $currentFeat) ;
		}
		elsif ($currentFeat->strand eq 1 and $oldFeat->strand eq -1) {
			$feature = &newSp_keepForward ($oldFeat, $currentFeat) ;
		}
		elsif ($currentFeat->strand eq -1 and $oldFeat->strand eq 1) {
			$feature = &newSp_keepForward ($oldFeat, $currentFeat) ;
		}
		return ($feature, $currentFeat) ;
	}
	elsif ( (($oldFeat->end-$oldFeat->start) / $oldFeat->{info}->{targetLgth} *100) >= 90 ) { # mite ou feature complete
		if ($currentFeat->strand eq $oldFeat->strand ) {
			$feature = &newSt_sameStrand ($currentFeat, $oldFeat) ;
		}
		elsif ($currentFeat->strand eq 1 and $oldFeat->strand eq -1) {
			$feature = &newSt_keepReverse ($currentFeat, $oldFeat) ;
		}
		elsif ($currentFeat->strand eq -1 and $oldFeat->strand eq 1) {
			$feature = &newSt_keepReverse ($currentFeat, $oldFeat) ;
		}
		return ($oldFeat, $feature) ;
	}
	elsif ((($currentFeat->end-$currentFeat->start) /  $currentFeat->{info}->{targetLgth} *100) >= 90 ) { # mite ou feature complete
		if ($currentFeat->strand eq $oldFeat->strand ) {
			$feature = &newSp_sameStrand ($oldFeat, $currentFeat) ;
		}
		elsif ($currentFeat->strand eq 1 and $oldFeat->strand eq -1) {
			$feature = &newSp_keepForward ($oldFeat, $currentFeat) ;
		}
		elsif ($currentFeat->strand eq -1 and $oldFeat->strand eq 1) {
			$feature = &newSp_keepForward ($oldFeat, $currentFeat) ;
		}
		return ($feature, $currentFeat) ;
	}
	elsif ($oldFeat->{info}->{targetEnd} + 50 > $oldFeat->{info}->{targetLgth} and  $currentFeat->{info}->{targetStart} < 50 ) { 
		goto (DOWN) ;
	}
	elsif (  $currentFeat->{info}->{targetEnd} + 50 >  $currentFeat->{info}->{targetLgth} and $oldFeat->{info}->{targetStart} < 50 ) { 
		goto (DOWN) ;
	}
	elsif ($oldFeat->{info}->{targetStart} < 50 ) {
		if ($currentFeat->strand eq $oldFeat->strand ) {
			$feature = &newSt_sameStrand ($currentFeat, $oldFeat) ;
		}
		elsif ($currentFeat->strand eq 1 and $oldFeat->strand eq -1) {
			$feature = &newSt_keepReverse ($currentFeat, $oldFeat) ;
		}
		elsif ($currentFeat->strand eq -1 and $oldFeat->strand eq 1) {
			$feature = &newSt_keepReverse ($currentFeat, $oldFeat) ;
		}
		return ($oldFeat, $feature) ;
	}
	elsif ($oldFeat->{info}->{targetEnd} + 50 > $oldFeat->{info}->{targetLgth} ) { # concerver une borne correct
		if ($currentFeat->strand eq $oldFeat->strand ) {
			$feature = &newSt_sameStrand ($currentFeat, $oldFeat) ;
		}
		elsif ($currentFeat->strand eq 1 and $oldFeat->strand eq -1) {
			$feature = &newSt_keepReverse ($currentFeat, $oldFeat) ;
		}
		elsif ($currentFeat->strand eq -1 and $oldFeat->strand eq 1) {
			$feature = &newSt_keepReverse ($currentFeat, $oldFeat) ;
		}
		return ($oldFeat, $feature) ;
	}

	elsif ( $currentFeat->{info}->{targetStart} < 50 ) {
		if ($currentFeat->strand eq $oldFeat->strand ) {
			$feature = &newSp_sameStrand ($oldFeat, $currentFeat) ;
		}
		elsif ($currentFeat->strand eq 1 and $oldFeat->strand eq -1) {
			$feature = &newSp_keepForward ($oldFeat, $currentFeat) ;
		}
		elsif ($currentFeat->strand eq -1 and $oldFeat->strand eq 1) {
			$feature = &newSp_keepForward ($oldFeat, $currentFeat) ;
		}
		return ($feature, $currentFeat) ;
	}

	elsif ( $currentFeat->{info}->{targetEnd} + 50 >  $currentFeat->{info}->{targetLgth} ) {
		if ($currentFeat->strand eq $oldFeat->strand ) {
			$feature = &newSp_sameStrand ($oldFeat, $currentFeat) ;
		}
		elsif ($currentFeat->strand eq 1 and $oldFeat->strand eq -1) {
			$feature = &newSp_keepForward ($oldFeat, $currentFeat) ;
		}
		elsif ($currentFeat->strand eq -1 and $oldFeat->strand eq 1) {
			$feature = &newSp_keepForward ($oldFeat, $currentFeat) ;
		}
		return ($feature, $currentFeat) ;
	}
	else {
		DOWN:if ($oldFeat->length >= $currentFeat->length) {
			if ($currentFeat->strand eq $oldFeat->strand ) {
				$feature = &newSt_sameStrand ($currentFeat, $oldFeat) ;
			}
			elsif ($currentFeat->strand eq 1 and $oldFeat->strand eq -1) {
				$feature = &newSt_keepReverse ($currentFeat, $oldFeat) ;
			}
			elsif ($currentFeat->strand eq -1 and $oldFeat->strand eq 1) {
				$feature = &newSt_keepReverse ($currentFeat, $oldFeat) ;
			}
			return ($oldFeat, $feature) ;
		}

		else {
			if ($currentFeat->strand eq $oldFeat->strand ) {
				$feature = &newSp_sameStrand ($oldFeat, $currentFeat) ;
			}
			elsif ($currentFeat->strand eq 1 and $oldFeat->strand eq -1) {
				$feature = &newSp_keepForward ($oldFeat, $currentFeat) ;
			}
			elsif ($currentFeat->strand eq -1 and $oldFeat->strand eq 1) {
				$feature = &newSp_keepForward ($oldFeat, $currentFeat) ;
			}
			return ($feature, $currentFeat) ;
		}
	}
}

sub newSt_sameStrand {
	my ($recalFeat, $keepFeat) = @_ ;
	my $newStart = $recalFeat->start + ($keepFeat->end - $recalFeat->start ) + 1 ;
	my $feature ;
	if ($recalFeat->strand eq 1) {
		my $newConsStart = $recalFeat->{info}->{targetStart} + ($keepFeat->end - $recalFeat->start ) + 1 ;
		if ( $newConsStart < 0 ) {$newConsStart = 1 ; }
		$feature = &createfeature ($recalFeat, $newStart, $recalFeat->end, $newConsStart , $recalFeat->{info}->{targetEnd}) ;
	}
	elsif ($recalFeat->strand eq -1) {
		my $newConsEnd = $recalFeat->{info}->{targetEnd} - ($keepFeat->end - $recalFeat->start ) + 1 ;
		if ( $newConsEnd < 0 ) {$newConsEnd = 1 ; }
		$feature = &createfeature ($recalFeat, $newStart, $recalFeat->end, $recalFeat->{info}->{targetStart}, $newConsEnd) ;
	}
	&checkFeat ($feature, "newSt_sameStrand") ;
	return ($feature) ;
}

sub newSp_sameStrand {
	my ($recalFeat, $keepFeat) = @_ ;
	my $newEnd = $recalFeat->end - ($recalFeat->end - $keepFeat->start ) - 1 ;
	my $feature ;
	if ($recalFeat->strand eq 1) {
		my $newConsEnd = $recalFeat->{info}->{targetEnd} - ($recalFeat->end - $keepFeat->start ) - 1 ;
		if ( $newConsEnd < 0 ) {$newConsEnd = 1 ; }
		$feature = &createfeature ($recalFeat, $recalFeat->start , $newEnd, $recalFeat->{info}->{targetStart}, $newConsEnd) ;
	}
	elsif ($recalFeat->strand eq -1) {
		my $newConsStart = $recalFeat->{info}->{targetStart} + ($recalFeat->end - $keepFeat->start ) + 1 ;
		if ( $newConsStart < 0 ) {$newConsStart = 1 ; }
		$feature = &createfeature ($recalFeat, $recalFeat->start , $newEnd, $newConsStart , $recalFeat->{info}->{targetEnd}) ;
	}
	&checkFeat ($feature, "newSp_sameStrand") ;
	return ($feature) ;
}

sub newSt_keepReverse {
	my ($recalFeat, $keepFeat) = @_ ;
	my $newStart = $recalFeat->start + ($keepFeat->end - $recalFeat->start ) + 1 ;
	my $feature ;
	if ($recalFeat->strand eq 1) {
		my $newConsStart = $recalFeat->{info}->{targetStart} + ($keepFeat->end - $recalFeat->start ) + 1 ;
		if ( $newConsStart < 0 ) {$newConsStart = 1 ; }
		$feature = &createfeature ($recalFeat, $newStart, $recalFeat->end, $newConsStart , $recalFeat->{info}->{targetEnd}) ;
	}
	elsif ($recalFeat->strand eq -1) {
		my $newConsEnd = $recalFeat->{info}->{targetEnd} - ($keepFeat->end - $recalFeat->start ) + 1 ;
		if ( $newConsEnd < 0 ) {$newConsEnd = 1 ; }
		$feature = &createfeature ($recalFeat, $newStart, $recalFeat->end, $recalFeat->{info}->{targetStart}, $newConsEnd) ;
	}
	&checkFeat ($feature, "newSt_keepReverse") ;
	return ($feature) ;
}

sub newSp_keepForward {
	my ($recalFeat, $keepFeat) = @_ ;
	my $newEnd = $recalFeat->end - ($recalFeat->end - $keepFeat->start ) - 1 ;
	my $feature ;
	if ($recalFeat->strand eq 1) {
		my $newConsEnd = $recalFeat->{info}->{targetEnd} - ($recalFeat->end - $keepFeat->start ) - 1 ;
		if ( $newConsEnd < 0 ) {$newConsEnd = 1 ; }
		$feature = &createfeature ($recalFeat, $recalFeat->start , $newEnd, $recalFeat->{info}->{targetStart}, $newConsEnd) ;
	}
	elsif ($recalFeat->strand eq -1) {
		my $newConsStart = $recalFeat->{info}->{targetStart} + ($recalFeat->end - $keepFeat->start ) + 1 ;
		if ( $newConsStart < 0 ) {$newConsStart = 1 ; }
		$feature = &createfeature ($recalFeat, $recalFeat->start , $newEnd, $newConsStart , $recalFeat->{info}->{targetEnd}) ;
	}
	&checkFeat ($feature, "newSp_keepForward") ;
	return ($feature) ;
}

sub checkFeat {
	my ($feature, $sub) = @_ ;
	if ($feature->{info}->{targetStart} < 0 or $feature->{info}->{targetEnd} < 0 or $feature->{info}->{targetEnd} < $feature->{info}->{targetStart}) {
		if ($verbosity ==4) { 
			print STDERR ("WARNING : wrong targetStart and/or targetEnd in new Feature from $sub (due to short feature): \n"
			."start\tend\ttargetStart\ttargetEnd\n"
			.$feature->start."\t".$feature->end."\t".$feature->{info}->{targetStart}."\t". $feature->{info}->{targetEnd}."\n") ;
		}
	}
	return 1 ;
}


#===== SUB Step 2 =====
sub filterOutSmallFeat {
	my ($i, $objFeat, $cutoff) = @_ ;
	if ($objFeat->{$i}->length < $cutoff
	and ($objFeat->{$i}->{info}->{targetStart} > 50 and $objFeat->{$i}->{info}->{targetEnd} + 50 < $objFeat->{$i}->{info}->{targetLgth})
	) {
		return 1 ;
	}
	return 0 ;
}

sub filterOutWongMatch {
	my ($prev, $cur, $next) = @_ ;
	my ($feat, $test) = &resolveLTR ($prev, $cur, $next) ;
	if ($test == 1 ) {
		$cur = $feat ;
	}
	my $newFeat ;
	my $cutoffPercent = 0 ;
	if ($prev->{info}->{targetFamID} =~/DT[A-Z]_fam.*/) { $cutoffPercent = 20 ; } # pour CACTA
	else {$cutoffPercent = 10 ; }
	
	if ($cur->{info}->{targetStart} > 50 and $cur->{info}->{targetEnd} + 50 < $cur->{info}->{targetLgth}
	and (($cur->end-$cur->start) / $cur->{info}->{targetLgth} *100) < $cutoffPercent # petite feature, sans borne
	and $prev->{info}->{targetFamID} eq $next->{info}->{targetFamID} # encadrée par 2 features appartenant à la même famille
	and $prev->{info}->{targetFamID} ne $cur->{info}->{targetFamID} # different des 2 feat qui l'encadre
		) {
			if ($prev->strand eq 1 and $next->strand eq 1
			and $prev->{info}->{targetEnd} + 50 < $prev->{info}->{targetLgth} and $next->{info}->{targetStart} > 50 # curfeat pas encadrée par sp et st
			and ( $prev->{info}->{targetEnd} < $next->{info}->{targetEnd} # colinarity
			or ( $prev->{info}->{targetStart} < $next->{info}->{targetStart} and $prev->{info}->{targetEnd} > $next->{info}->{targetEnd} )) # include
			) {
				$newFeat = &mergefeature ($prev, $next, $prev->start, $next->end, $prev->{info}->{targetStart}, $next->{info}->{targetEnd} ) ;
				return ($newFeat) ;
			}
			elsif ($prev->strand eq -1 and $next->strand eq -1
			and $prev->{info}->{targetStart} > 50 and $next->{info}->{targetEnd} + 50 < $next->{info}->{targetLgth}  # curfeat pas encadrée par sp et st
			and ($prev->{info}->{targetStart} > $next->{info}->{targetStart} # colinarity
			or ($prev->{info}->{targetEnd} > $next->{info}->{targetEnd} and $prev->{info}->{targetStart} > $next->{info}->{targetStart} )) # include
			) {
				$newFeat = &mergefeature ($prev, $next, $prev->start, $next->end, $next->{info}->{targetStart}, $prev->{info}->{targetEnd} ) ;
				return ($newFeat) ;
			}
			else { return 0 ; } 
		}
	else { return 0 ; }
}

sub mergeColFeat {
	my ($prev, $cur, $next, $self) = @_ ;
	if ($next != 0 ) {
		my ($feat, $test) = &resolveLTR ($prev, $cur, $next, $self) ;
		if ($test == 1 ) {
			$cur = $feat ;
			if ($verbosity == 4) { print STDERR join ("\t", "LTRresolve", $cur->{info}->{targetVarID}, $cur->start, $cur->strand, $cur->end, $cur->{info}->{targetStart}, $cur->{info}->{targetEnd}), "\n" ; }
		}
	}
	my $return ;
	my $newFeat ;
	if ($cur->{info}->{targetFamID} eq $prev->{info}->{targetFamID} and $cur->strand eq $prev->strand) {
		if ($cur->strand eq 1 and $prev->strand eq 1) { # brin plus
			if ( defined ($self->{_classification}->{$prev->{info}->{copieID}}->{"LTR5"}->{end}) and defined ($self->{_classification}->{$cur->{info}->{copieID}}->{"LTR5"}->{end}) # soloLTR : prev LTR5' et cur LTR3' => Merge
			and $prev->{info}->{targetEnd}  < $self->{_classification}->{$prev->{info}->{copieID}}->{"LTR5"}->{end}
			and $cur->{info}->{targetStart} > $self->{_classification}->{$cur->{info}->{copieID}}->{"LTR3"}->{start} ) {
				if ($cur->start - $prev->end < 200) {
					$newFeat = &mergefeature ($prev, $cur, $prev->start, $cur->end, $prev->{info}->{targetStart}, $cur->{info}->{targetEnd} ) ;
					$newFeat->{soloLTR} = 1 ;
					$return->{newFeat} = $newFeat ;
				}
				else {$return->{cur} = $cur ; }
			}
			elsif ($prev->{info}->{targetEnd} < $cur->{info}->{targetEnd}  # colinearity of preidction => Merge
			and $prev->{info}->{targetEnd} + 50 < $prev->{info}->{targetLgth}
			and $cur->end - $prev->start < (($prev->{info}->{targetLgth} + $cur->{info}->{targetLgth})/2) * 1.5
			and $cur->{info}->{targetStart} > 50 ) {
				$newFeat = &mergefeature ($prev, $cur, $prev->start, $cur->end, $prev->{info}->{targetStart}, $cur->{info}->{targetEnd} ) ;
				$return->{newFeat} = $newFeat ;
			}
			elsif (
			$prev->{info}->{targetEnd} > $cur->{info}->{targetEnd} and $prev->{info}->{targetStart} < $cur->{info}->{targetStart} # cur inclue in prev => Merge
			and $cur->end - $prev->start < (($prev->{info}->{targetLgth} + $cur->{info}->{targetLgth})/2) * 1.5
			and $prev->{info}->{targetEnd} + 50 < $prev->{info}->{targetLgth} ) {
				$newFeat = &mergefeature ($prev, $cur, $prev->start, $cur->end, $prev->{info}->{targetStart}, $prev->{info}->{targetEnd} ) ;
				$return->{newFeat} = $newFeat ;
			}
			else { $return->{cur} = $cur ; }
		}
		elsif ($cur->strand eq -1 and $prev->strand eq -1) { # brin moins
			if ( defined ($self->{_classification}->{$prev->{info}->{copieID}}->{"LTR3"}->{start}) and defined ($self->{_classification}->{$cur->{info}->{copieID}}->{"LTR5"}->{end}) # soloLTR
			and $prev->{info}->{targetStart}  > $self->{_classification}->{$prev->{info}->{copieID}}->{"LTR3"}->{start}
			and $cur->{info}->{targetEnd} < $self->{_classification}->{$cur->{info}->{copieID}}->{"LTR5"}->{end} ) { 
				if ($cur->end - $prev->start < 200) {
					$newFeat = &mergefeature ($prev, $cur, $prev->start, $cur->end, $cur->{info}->{targetStart}, $prev->{info}->{targetEnd} ) ;
					$newFeat->{soloLTR} = 1 ;
					$return->{newFeat} = $newFeat ;
				}
				else { $return->{cur} = $cur ; }
				}

			elsif ($prev->{info}->{targetStart} > $cur->{info}->{targetStart} # collinear
			and $prev->{info}->{targetStart} > 50
			and $cur->end - $prev->start < (($prev->{info}->{targetLgth} + $cur->{info}->{targetLgth})/2) * 1.5
			and $cur->{info}->{targetEnd} + 50 < $cur->{info}->{targetLgth} ) {
				$newFeat = &mergefeature ($prev, $cur, $prev->start, $cur->end, $cur->{info}->{targetStart}, $prev->{info}->{targetEnd} ) ;
				$return->{newFeat} = $newFeat ;
			}
			elsif (
			$prev->{info}->{targetStart} < $cur->{info}->{targetStart} and $prev->{info}->{targetEnd} > $cur->{info}->{targetEnd} # inclue
			and $cur->end - $prev->start < (($prev->{info}->{targetLgth} + $cur->{info}->{targetLgth})/2) * 1.5
			and $prev->{info}->{targetStart} > 50 ) {
				$newFeat = &mergefeature ($prev, $cur, $prev->start, $cur->end, $prev->{info}->{targetStart}, $prev->{info}->{targetEnd} ) ;
				$return->{newFeat} = $newFeat ;
			}
			else { $return->{cur} = $cur ; }
		}
	}
	else { $return->{cur} = $cur ; }
	return ($return) ;
}

sub resolveLTR {
	my ($prev, $cur, $next, $self) = @_ ;
	my $newFeat ;
	my $test = 0 ;
	if (defined($self->{_classification}->{$cur->{info}->{copieID}}->{"LTR5"}->{end})
	and $cur->{info}->{targetEnd} <= $self->{_classification}->{$cur->{info}->{copieID}}->{"LTR5"}->{end} + 500 ) { # séléction des feature inclue dans le LTR 5' + 500 pb
		if ($cur->strand eq 1 and $prev->strand eq 1 ) {
			($newFeat, $test) = &wrongLTR5 ($prev, $cur, $next, $self) ;
			if ($test == 0 ) { ($newFeat, $test) = &wrongLTR5Inside ($prev, $cur, $next, $self) ; }
		}
		elsif ($cur->strand eq -1 and $next->strand eq -1 ) {
			($newFeat, $test) = &wrongLTR5 ($next, $cur, $prev, $self) ;
			if ($test == 0 ) { ($newFeat, $test) = &wrongLTR5Inside ($prev, $cur, $next, $self) ; }
		}
	}
	elsif (defined($self->{_classification}->{$cur->{info}->{copieID}}->{"LTR3"}->{end})
	and $cur->{info}->{targetStart} >= $self->{_classification}->{$cur->{info}->{copieID}}->{"LTR3"}->{start} - 500 ) { # séléction des feature inclue dans le LTR 3' + 500 pb
		if ($cur->strand eq 1 and $next->strand eq 1) {
			($newFeat, $test) = &wrongLTR3 ($prev, $cur, $next, $self) ;
			if ($test == 0 ) { ($newFeat, $test) = &wrongLTR3Inside ($prev, $cur, $next, $self) ; }
		}
		elsif ($cur->strand eq -1 and $prev->strand eq -1) {
			($newFeat, $test) = &wrongLTR3 ($next, $cur, $prev, $self) ;
			if ($test == 0 ) { ($newFeat, $test) = &wrongLTR3Inside ($prev, $cur, $next, $self) ; }
		}
	}
	return ($newFeat, $test) ;
}

sub wrongLTR5Inside { # preidction of a LTR5 are between two LTR3 => LTR3 LTR5 LTR3
	my ($prev, $wrong, $next, $self) = @_ ;
	my $newFeat ;
	my $test = 0 ;
	if ( ($wrong->{info}->{targetFamID} eq $prev->{info}->{targetFamID} and $prev->{info}->{targetFamID} eq $next->{info}->{targetFamID})
	and ($wrong->strand eq $prev->strand and $prev->strand eq $next->strand)
	 ) {
		if ( $wrong->strand eq 1
		and defined($self->{_classification}->{$next->{info}->{copieID}}->{"LTR3"}->{end})
		and $next->{info}->{targetStart} > $self->{_classification}->{$next->{info}->{copieID}}->{"LTR3"}->{start}
		and ( abs(($next->{info}->{targetStart} - $prev->{info}->{targetEnd}) - ($next->start - $prev->end) ) < 500)
		and ($next->end - $prev->start < $prev->{info}->{targetLgth} + 1000)
		){
			$newFeat = &createfeature ($prev, $wrong->start, $wrong->end, ($prev->{info}->{targetEnd} + 1) , ($next->{info}->{targetStart} - 1) ) ;
			$test = 1 ;
		}
		elsif ( $wrong->strand eq -1
		and defined($self->{_classification}->{$prev->{info}->{copieID}}->{"LTR3"}->{end})
		and $prev->{info}->{targetStart} > $self->{_classification}->{$prev->{info}->{copieID}}->{"LTR3"}->{start}
		and ( abs(abs($prev->{info}->{targetStart} - $next->{info}->{targetEnd}) - ($next->start - $prev->end)) < 500)
		and ($prev->start - $next->end < $prev->{info}->{targetLgth} + 1000)
		){
			$newFeat = &createfeature ($prev, $wrong->start, $wrong->end, ($next->{info}->{targetEnd} + 1), ($prev->{info}->{targetStart} - 1) ) ;
			$test = 1 ;
		}
	}
	return ($newFeat, $test) ;
}

sub wrongLTR3Inside { # match of a LTR3' are between two LTR5' => LTR5 LTR3 LTR5  
	my ($prev, $wrong, $next, $self) = @_ ;
	my $newFeat ;
	my $test = 0 ;
	if ( ($wrong->{info}->{targetFamID} eq $prev->{info}->{targetFamID} and $prev->{info}->{targetFamID} eq $next->{info}->{targetFamID})
	and ( $wrong->strand eq $prev->strand and $prev->strand eq $next->strand )
	 ) {
		if ( $wrong->strand eq 1
		and defined($self->{_classification}->{$prev->{info}->{copieID}}->{"LTR5"}->{end})
		and $prev->{info}->{targetEnd} < $self->{_classification}->{$prev->{info}->{copieID}}->{"LTR5"}->{end}
		and ( abs(($next->{info}->{targetStart} - $prev->{info}->{targetEnd}) - ($next->start - $prev->end) ) < 500)
		and ($next->end - $prev->start < $prev->{info}->{targetLgth} + 1000)
		){
			$newFeat = &createfeature ($next, $wrong->start, $wrong->end, ($prev->{info}->{targetEnd} + 1) , ($next->{info}->{targetStart} - 1) ) ;
			$test = 1 ;
		}
		elsif ( $wrong->strand eq -1 
		and defined($self->{_classification}->{$prev->{info}->{copieID}}->{"LTR5"}->{end})
		and $next->{info}->{targetEnd} < $self->{_classification}->{$prev->{info}->{copieID}}->{"LTR5"}->{end} # prédiction n'est pas à la fin du LTR 5
		and ( abs(($prev->{info}->{targetStart} - $next->{info}->{targetEnd}) - ($next->start- $prev->end) < 500))
		and ($next->start - $prev->end < $prev->{info}->{targetLgth} + 1000)
		){
				$newFeat = &createfeature ($next, $wrong->start, $wrong->end, ($next->{info}->{targetEnd} + 1), ($prev->{info}->{targetStart} - 1) ) ;
				$test = 1 ;
			}
	}
	return ($newFeat, $test) ;
}

sub wrongLTR5 { # wrong LTR5 at the end of the RT-LTR : LTR5-TE-LTR5
	my ($true, $wrong, $rand, $self) = @_ ;
	my $newFeat ;
	my $test = 0 ;
	if (defined($self->{_classification}->{$true->{info}->{copieID}}->{"LTR3"}->{end})) {
		if ( ($wrong->{info}->{targetFamID} eq $true->{info}->{targetFamID} and $true->{info}->{targetFamID} ne $rand->{info}->{targetFamID})
		or ( $wrong->strand eq $true->strand and $true->strand ne $rand->strand and $wrong->{info}->{targetFamID} eq $true->{info}->{targetFamID} )
		 ) {
			if ( $wrong->strand eq 1
			and (abs ( abs ($wrong->end - $true->end) - ($true->{info}->{targetLgth} - $true->{info}->{targetEnd} ) ) < 500 
				 or abs ( abs ($wrong->{info}->{targetEnd} - $wrong->{info}->{targetStart}) - ($true->{info}->{targetLgth} - $true->{info}->{targetEnd} ) ) < 500)
			){
				$newFeat = &createfeature ($true, $wrong->start, $wrong->end, $self->{_classification}->{$true->{info}->{copieID}}->{"LTR3"}->{start}, $self->{_classification}->{$true->{info}->{copieID}}->{"LTR3"}->{end}) ;
				$test = 1 ;
			}
			elsif ( $wrong->strand eq -1
			and (abs ( abs ($wrong->start - $true->start) - ($true->{info}->{targetLgth} - $true->{info}->{targetEnd} ) ) < 500 
				or abs ( abs ($wrong->{info}->{targetEnd} - $wrong->{info}->{targetStart}) - ($true->{info}->{targetLgth} - $true->{info}->{targetEnd} ) ) < 500)
			){
				$newFeat = &createfeature ($true, $wrong->start, $wrong->end, $self->{_classification}->{$true->{info}->{copieID}}->{"LTR3"}->{start}, $self->{_classification}->{$true->{info}->{copieID}}->{"LTR3"}->{end}) ;
				$test = 1 ;
			}
		}
	}
	return ($newFeat, $test) ;
}

sub wrongLTR3 { # wrong LTR3 at the beginning of the RT-LTR : LTR3-TE-LTR3
	my ($rand, $wrong, $true, $self) = @_ ;
	my $test = 0 ;
	my $newFeat ;
	if (defined($self->{_classification}->{$true->{info}->{copieID}}->{"LTR5"}->{end})) {
		if ( ($wrong->{info}->{targetFamID} eq $true->{info}->{targetFamID} and $true->{info}->{targetFamID} ne $rand->{info}->{targetFamID})
		or ($wrong->strand eq $true->strand and $true->strand ne $rand->strand and $wrong->{info}->{targetFamID} eq $true->{info}->{targetFamID} )
		){
			if ( $wrong->strand eq 1
			and $true->{info}->{targetStart} > 50 
			and (abs ( ($true->start - $wrong->start) - $true->{info}->{targetStart} ) < 500 or abs ( ($wrong->{info}->{targetEnd} - $wrong->{info}->{targetStart}) - $true->{info}->{targetStart} ) < 500) 
			){
				$newFeat = &createfeature ($true, $wrong->start, $wrong->end, $self->{_classification}->{$true->{info}->{copieID}}->{"LTR5"}->{start}, $self->{_classification}->{$true->{info}->{copieID}}->{"LTR5"}->{end}) ;
				$test = 1 ;
			}
			elsif ($wrong->strand eq -1
			and $true->{info}->{targetStart} > 50 
			and (abs ( $true->{info}->{targetStart} - ($wrong->end - $true->end)) < 500 or abs ( $true->{info}->{targetStart} - ($wrong->{info}->{targetEnd} - $wrong->{info}->{targetStart})) < 500)
			){
				$newFeat = &createfeature ($true, $wrong->start, $wrong->end, $self->{_classification}->{$true->{info}->{copieID}}->{"LTR5"}->{start}, $self->{_classification}->{$true->{info}->{copieID}}->{"LTR5"}->{end}) ;
				$test = 1 ;
			}
		}
	}
	return ($newFeat, $test) ;
}

sub mergefeature {
	my ($obj_1, $obj_2, $start, $end, $consStart, $consEnd) = @_ ;
	my $seqID = $obj_1->seq_id ;
	my $pritag = "repeat_region" ;
	my $location = $obj_1->location ;
	my $strand = $location->strand ;
	my $id = ($obj_1 ->each_tag_value("id"))[0] ;
	my $loc = Bio::Location::Split->new ;
	my $feature = Bio::SeqFeature::Generic -> new ( -seq_id      => $seqID,
													-source_tag  => $source_tag,
													-primary_tag => $pritag,
													-start       => $start,
													-end         => $end,
													-strand      => $strand
													 );
													 
	foreach my $variant (keys %{$obj_1->{compo}}) {
		$feature->{compo}->{$variant}->{lgth} = $obj_1->{compo}->{$variant}->{lgth} ;
		$feature->{compo}->{$variant}->{conslgth} = $obj_1->{compo}->{$variant}->{conslgth} ;
	}
	$feature->{compo}->{$obj_2->{info}->{targetVarID}}->{lgth} += $obj_2->length ;
	$feature->{compo}->{$obj_2->{info}->{targetVarID}}->{conslgth} = $obj_2->{info}->{targetLgth} ;
	my $testEnd = 0 ;
	if ($obj_2->strand eq "1" and $obj_2->{info}->{targetEnd} + 50 > $obj_2->{info}->{targetLgth}) { # case where the feature obj_2 is at the extemity of the reference sequence
		$testEnd = 1  
	}
	elsif ($obj_1->strand eq "-1" and $obj_1->{info}->{targetEnd} + 50 > $obj_1->{info}->{targetLgth}) {
		$testEnd = 1  
	}
	my $var = 0 ;
	my $mergeK = 0 ;
	foreach my $allvariant (keys %{$feature->{compo}}) {
		$mergeK++ ;
		if ($mergeK == 1 ) { $var = $allvariant ; next ; }
		if ($feature->{compo}->{$var}->{lgth} < $feature->{compo}->{$allvariant}->{lgth}) { $var = $allvariant ; }
	}
	my $target ;
	if ($testEnd == 0) { $target = $var." ".$feature->{compo}->{$var}->{conslgth}."bp ".$consStart."..".$consEnd ; $feature->{info}->{targetEnd} = $consEnd ; }
	if ($testEnd == 1) { $target = $var." ".$feature->{compo}->{$var}->{conslgth}."bp ".$consStart."..".$feature->{compo}->{$var}->{conslgth} ; $feature->{info}->{targetEnd} = $feature->{compo}->{$var}->{conslgth} ; }
	$feature ->add_tag_value ( "post", $target ) ;
	$feature ->add_tag_value ( "id", $id ) ;
	
	$feature->{info}->{targetVarID} = $var ;
	if ($var=~/^([A-Z]*_fam[n|c][0-9]*)\.[0-9]+$/) { $feature->{info}->{targetFamID} = $1 ; }
	elsif ($var=~/^([A-Z]*_fam[n|c][0-9]+)$/) {	$feature->{info}->{targetFamID} = $1 ; }
	$feature->{info}->{targetStart} = $consStart ;
	$feature->{info}->{targetLgth} = $feature->{compo}->{$var}->{conslgth} ;
	if ($obj_1->length >= $obj_2->length) {$feature->{info}->{copieID} = $obj_1->{info}->{copieID} ; }
	elsif ($obj_1->length < $obj_2->length) {$feature->{info}->{copieID} = $obj_2->{info}->{copieID} ; }
	return ($feature) ;
}

sub createfeature {
	my ($obj_feat, $start, $end, $consStart, $consEnd) = @_ ;
	my $seqID = $obj_feat->seq_id ;
	my $pritag = "repeat_region" ;
	my $location = $obj_feat->location ;
	my $strand = $location->strand ;
	my $ID = ($obj_feat ->each_tag_value("id"))[0] ;
	my $loc = Bio::Location::Split->new ;
	my $feature = Bio::SeqFeature::Generic -> new ( -seq_id      => $seqID,
													-source_tag  => $source_tag,
													-primary_tag => $pritag,
													-start       => $start,
													-end         => $end,
													-strand      => $strand
													 ); 
	$feature ->add_tag_value ( "id", $ID ) ;
	my $target = $obj_feat->{info}->{targetVarID}." ".$obj_feat->{info}->{targetLgth}."bp ".$consStart."..".$consEnd ;
	$feature ->add_tag_value ( "post", $target ) ;
	$feature->{compo}->{$obj_feat->{info}->{targetVarID}}->{lgth} = $feature->length ;
	$feature->{compo}->{$obj_feat->{info}->{targetVarID}}->{conslgth} = $obj_feat->{info}->{targetLgth} ;
	$feature->{info}->{targetVarID} = $obj_feat->{info}->{targetVarID};
	$feature->{info}->{targetFamID} = $obj_feat->{info}->{targetFamID};
	$feature->{info}->{targetStart} = $consStart ;
	$feature->{info}->{targetEnd} = $consEnd ;
	$feature->{info}->{targetLgth} = $obj_feat->{info}->{targetLgth};
	$feature->{info}->{copieID} = $obj_feat->{info}->{copieID} ;
	return ($feature) ;
}

sub createParentFeature {
	my ($obj_feat1, $obj_feat2) = @_ ;
	my $seqID = ${$obj_feat1}->seq_id ;
	my $pritag = "match_part" ;
	my $location1 = ${$obj_feat1}->location ;
	my $location2 = ${$obj_feat2}->location ;
	my $strand = $location1->strand ;
	my $start  = $location1->start ;
	my $end    = $location2->end ;
	my $note ;
	my $feature = Bio::SeqFeature::Generic -> new ( -seq_id      => $seqID,
													-source_tag  => $source_tag,
													-primary_tag => $pritag,
													-start       => $start,
													-end         => $end,
													-strand      => $strand
													 );
	# parent
	my $parent ;
	if (${$obj_feat1}->primary_tag eq "match_part" ) { $parent = (${$obj_feat1}->each_tag_value("id"))[0] ; }
	elsif (${$obj_feat2}->primary_tag eq "match_part" ) { $parent = (${$obj_feat2}->each_tag_value("id"))[0] ; }
	else { 
		my $matchId = (${$obj_feat1}->each_tag_value("id"))[0] ; 
		$kp++ ;
		$parent = "mp".$kp ;
	}
	
	# range
	if (${$obj_feat1}->{info}->{targetFamID} ne "exon" and ${$obj_feat2}->{info}->{targetFamID} ne "exon") {
		if (${$obj_feat1}->primary_tag eq "match_part" ) { 
			foreach my $targStart (keys %{${$obj_feat1}->{range}}){
				$feature->{range}->{$targStart} = ${$obj_feat1}->{range}->{$targStart}
			}
		}
		else { $feature->{range}->{${$obj_feat1}->{info}->{targetStart}} ="${$obj_feat1}->{info}->{targetStart}..${$obj_feat1}->{info}->{targetEnd}" ; }
		if (${$obj_feat2}->primary_tag eq "match_part" ) { 
			foreach my $targStart (keys %{${$obj_feat2}->{range}}){
				$feature->{range}->{$targStart} = ${$obj_feat2}->{range}->{$targStart}
			}
		}
		else { $feature->{range}->{${$obj_feat2}->{info}->{targetStart}} ="${$obj_feat2}->{info}->{targetStart}..${$obj_feat2}->{info}->{targetEnd}" ; }
	}
	
	
	# composition
	foreach my $variant (keys %{${$obj_feat1}->{compo}}) {
		$feature->{compo}->{$variant}->{lgth} = ${$obj_feat1}->{compo}->{$variant}->{lgth} ;
		$feature->{compo}->{$variant}->{conslgth} = ${$obj_feat1}->{compo}->{$variant}->{conslgth} ;
	}
	foreach my $variant (keys %{${$obj_feat2}->{compo}}) {
		$feature->{compo}->{$variant}->{lgth} += ${$obj_feat2}->{compo}->{$variant}->{lgth} ;
		$feature->{compo}->{$variant}->{conslgth} = ${$obj_feat2}->{compo}->{$variant}->{conslgth} ;
	}
	my $testEnd = 0 ;
	if (${$obj_feat2}->strand eq "1" and ${$obj_feat2}->{info}->{targetEnd} + 50 > ${$obj_feat2}->{info}->{targetLgth}) { # cas ou la feature obj_2 et à l'extrémité du consensus
		$testEnd = 1  
	}
	elsif (${$obj_feat1}->strand eq "-1" and ${$obj_feat1}->{info}->{targetEnd} + 50 > ${$obj_feat1}->{info}->{targetLgth}) {
		$testEnd = 1  
	}
	my $var = 0 ;
	my $mergeK = 0 ;
	foreach my $allvariant (keys %{$feature->{compo}}) {
		$mergeK++ ;
		if ($mergeK == 1 ) { $var = $allvariant ; next ; }
		if ($feature->{compo}->{$var}->{lgth} < $feature->{compo}->{$allvariant}->{lgth}) { $var = $allvariant ; }
	}
	my $target ;
	if ($strand eq 1) {
		if ($testEnd == 0) {
			$target = $var." ".$feature->{compo}->{$var}->{conslgth}."bp ".${$obj_feat1}->{info}->{targetStart}."..".${$obj_feat2}->{info}->{targetEnd} ; 
			$feature->{info}->{targetStart} = ${$obj_feat1}->{info}->{targetStart} ;
			$feature->{info}->{targetEnd} = ${$obj_feat2}->{info}->{targetEnd} ;
		}
		if ($testEnd == 1) {
			$target = $var." ".$feature->{compo}->{$var}->{conslgth}."bp ".${$obj_feat1}->{info}->{targetStart}."..".$feature->{compo}->{$var}->{conslgth} ; 
			$feature->{info}->{targetStart} = ${$obj_feat1}->{info}->{targetStart} ;
			$feature->{info}->{targetEnd} = $feature->{compo}->{$var}->{conslgth} ;
		}
	}
	else {
		if ($testEnd == 0) {
			$target = $var." ".$feature->{compo}->{$var}->{conslgth}."bp ".${$obj_feat2}->{info}->{targetStart}."..".${$obj_feat1}->{info}->{targetEnd} ; 
			$feature->{info}->{targetStart} = ${$obj_feat2}->{info}->{targetStart} ;
			$feature->{info}->{targetEnd} = ${$obj_feat1}->{info}->{targetEnd} ;
		}
		if ($testEnd == 1) {
			$target = $var." ".$feature->{compo}->{$var}->{conslgth}."bp ".${$obj_feat2}->{info}->{targetStart}."..".$feature->{compo}->{$var}->{conslgth} ;
			$feature->{info}->{targetStart} = ${$obj_feat2}->{info}->{targetStart} ;
			$feature->{info}->{targetEnd} = $feature->{compo}->{$var}->{conslgth} ;
		}
	}
	${$obj_feat1} ->add_tag_value ( "parent", $parent ) ;
	${$obj_feat2} ->add_tag_value ( "parent", $parent ) ;
	$feature ->add_tag_value ( "id", $parent ) ;
	$feature ->add_tag_value ( "post", $target ) ;
	$feature->{info}->{targetVarID} = $var ;
	if ($var=~/^([A-Z]*_fam[n|c][0-9]*)\.[0-9]+$/) { $feature->{info}->{targetFamID} = $1 ; }
	elsif ($var=~/^([A-Z]*_fam[n|c][0-9]+)$/) {$feature->{info}->{targetFamID} = $1 ; }
	$feature->{info}->{targetLgth} = $feature->{compo}->{$var}->{conslgth} ;
	if (${$obj_feat1}->length >= ${$obj_feat2}->length) {$feature->{info}->{copieID} = ${$obj_feat1}->{info}->{copieID} ; }
	elsif (${$obj_feat1}->length < ${$obj_feat2}->length) {$feature->{info}->{copieID} = ${$obj_feat2}->{info}->{copieID} ; }
	return ($feature) ;
}

sub overlaping {
	my $self = shift ;
	my $oldFeat ;
	my $kfeat = 0 ;
	for (my $i = 1 ; $i <= scalar (keys (%{$self->{feat}})) ; $i++ ) {
		if ($i == 1 ) { $oldFeat = $self->{feat}->{$i} ; next ; }
		$kfeat++ ;
		TOP:
		if (not defined ($oldFeat)) { print STDERR ( "*** old obj is undefined value\n" ) and die ; }
		if (not defined ($self->{feat}->{$i})) { print STDERR ( "*** current obj is undefined value" ) and die ; }
		if ( $oldFeat->end >= $self->{feat}->{$i}->start 
		and $oldFeat->end < $self->{feat}->{$i}->end and $oldFeat->start < $self->{feat}->{$i}->start ) { # is overlapping but not include
			my ($push, $tmp) = &resolveOverlap ($oldFeat, $self->{feat}->{$i} ) ; 
			$self->{noOverlapingFeat}->{$kfeat} = $push ;
			$oldFeat = $tmp ;
			if ($i == scalar (keys (%{$self->{feat}})) ) { $self->{noOverlapingFeat}->{$kfeat+1} = $tmp ; }
		}
		elsif ( $oldFeat->start >= $self->{feat}->{$i}->start and $oldFeat->end <= $self->{feat}->{$i}->end ) { # old feat is include in next feat
			$kfeat = $kfeat-1 ;
			$oldFeat = $self->{noOverlapingFeat}->{$kfeat} ; 
			if ($kfeat == 0 and $i < scalar (keys (%{$self->{feat}})) ) { $kfeat=$kfeat+1 ; $oldFeat = $self->{feat}->{$i} ; $i=$i+1 ; }
			if ( $i == scalar (keys (%{$self->{feat}})) ) { $self->{noOverlapingFeat}->{$kfeat} = $self->{feat}->{$i} ; last }
			goto (TOP) ;
		}
		elsif ( $oldFeat->start <= $self->{feat}->{$i}->start and $oldFeat->end >= $self->{feat}->{$i}->end ) { # next feat is include in old feat 
			$kfeat = $kfeat-1 ;
		}
		elsif ($self->{feat}->{$i}->start < $oldFeat->start and $oldFeat->end >= $self->{feat}->{$i}->end ) { # old start > cur start, infertion de oldfeat et curfeat
			$kfeat = $kfeat-1 ;
			$oldFeat = $self->{noOverlapingFeat}->{$kfeat} ;
			$self->{feat}->{$i-1} = $self->{feat}->{$i} ;
			$self->{feat}->{$i} = $oldFeat ;
			$i = $i - 1 ;
			goto (TOP) ;
		}
		
		else { # no overlap between features
			$self->{noOverlapingFeat}->{$kfeat} = $oldFeat ;
			$oldFeat = $self->{feat}->{$i} ;
			if ($i == scalar (keys (%{$self->{feat}})) ) { $self->{noOverlapingFeat}->{$kfeat+1} = $oldFeat ; }
		}
	}
}

#===== SUB Step 3 =====

sub joinFeat {
	my ($i, $j, $mergeFeat, $threshold) = @_ ;
	my $test = 0 ;
	if ($mergeFeat->{$i+$j}->{info}->{targetFamID} eq $mergeFeat->{$i-1}->{info}->{targetFamID}
	and $mergeFeat->{$i-1}->{info}->{targetFamID} ne "gap"
	and (($mergeFeat->{$i-1}->strand eq 1 and $mergeFeat->{$i+$j}->strand eq 1 and $mergeFeat->{$i-1}->{info}->{targetEnd} + 10 < $mergeFeat->{$i-1}->{info}->{targetLgth} and $mergeFeat->{$i+$j}->{info}->{targetStart} > 10)
	or ( $mergeFeat->{$i-1}->strand eq -1 and $mergeFeat->{$i+$j}->strand eq -1 and $mergeFeat->{$i+$j}->{info}->{targetEnd} + 10 < $mergeFeat->{$i+$j}->{info}->{targetLgth} and $mergeFeat->{$i-1}->{info}->{targetStart} > 10))
	and (( $mergeFeat->{$i-1}->strand eq 1 and $mergeFeat->{$i+$j}->strand eq 1 and abs ( $mergeFeat->{$i+$j}->{info}->{targetStart} - $mergeFeat->{$i-1}->{info}->{targetEnd}) < $threshold )
	or ( $mergeFeat->{$i-1}->strand eq -1 and $mergeFeat->{$i+$j}->strand eq -1 and abs ( $mergeFeat->{$i-1}->{info}->{targetStart} - $mergeFeat->{$i+$j}->{info}->{targetEnd}) < $threshold ))
	) {
		$test = 1 ;
	}
	return ($test) ;
}

#==================================================================================
sub help {
my $prog = basename($0) ;
print STDERR <<EOF ;
#### $prog ####
#
# CREATED:    2012-05-31
# LAST MODIF: $lastmodif
# AUTHOR:     Josquin Daron (INRA Clermont-Ferrand)
# VERSION:    $VERSION
#
# This script is used to post process a embl output of RepeatMasker
# 

USAGE:
       $prog -fasta <fasta_file> -LTR <position of LTR> -classi <classification> -gene <embl> <xm file>


       ### OPTIONS ###

       -h|--help:             print this help
       -LTR:                  tabulation file of the LTR annotations
       -classi:               corresponding file sequence id to family id
       -fasta:                fasta file of the annotated sequence
       -gene:                 embl file of the gene annotation
       -dir:                  directory name
       -v:                    verbosity (3,4)

EOF
exit(1) ;
}
#==================================================================================
