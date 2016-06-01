#!/usr/bin/perl
#program which determine the list and category of each pseudogene candidate using the output of the first module
#Psi-Fi module 2 - Lerat and Ochman (2004) Genome Res.
#to use with the blast+ suite where 
#to construct the database of sequences
#
#$makeblastdb -in inputfile -dbtype 'nucl' -title databasetitle -input_type 'type' - out database_name -parse_seqids
#
#to obtain tabular output format
#tblastn -db database_name -query input_file -out output_file -outfmt 6

use Bio::SeqIO;

print "************************ Psi_Fi_module2 **********************\n";

print "Enter the name of the file\n";
$file = <STDIN>;
chomp($file);
open (FICHIER, "<$file") || die "File not found ";
@candidate=<FICHIER>;
close(FICHIER);

print "Enter the name of the target genome file (genbank format)\n";
$fileb = <STDIN>;
chomp($fileb);

my $seqio_object = Bio::SeqIO->new(-format => 'genbank', -file => "$fileb");
my $seq_object = $seqio_object->next_seq;

print "Enter the name of the out file\n";
$outfile = <STDIN>;
chomp($outfile);
open (OUTFILE, ">$outfile");

print "Enter the name of the genome bank (blast database name)\n";
$bank = <STDIN>;
chomp($bank);

$evalue=1e-15;
$prop_target=0.8;
$prop_match_min=0.3;
$prop_match_max=0.8;

print "************************ Parameters **********************\n";
print "A\tmaximal e_value: $evalue\n";
print "B\tmaximal percentage of the query size: $prop_target\n";
print "C\tminimal percentage of the query size (match): $prop_match_min\n";
print "D\tmaximal percentage of the query size (match): $prop_match_max\n";
print "Select the parameter you want to change or type Y\n";
$response=<STDIN>;
chomp($response);

while($response ne "Y")
{
	if($response eq "A")
	{
		print "Enter new value for e_value\n";
		$evalue=<STDIN>;
		chomp($evalue);
	}
	elsif($response eq "B")	
	{
		print "Enter new value for maximal length of the query size\n";
		$prop_target=<STDIN>;
		chomp($prop_target);
	}
	elsif($response eq "C")
	{
		print "Enter new value for minimal length of the query size (match)\n";
		$prop_match_min=<STDIN>;
		chomp($prop_match_min);
	}
	elsif($response eq "D")
	{
		print "Enter new value for maximal length of the query size (match)\n";
		$prop_match_max=<STDIN>;
		chomp($prop_match_max);
	}

	else
	{
		print "character not valid\n";
	}
	print "************************ Parameters **********************\n";
	print "A\tmaximal e_value: $evalue\n";
	print "B\tmaximal percentage of the query size: $prop_target\n";
	print "C\tminimal percentage of the query size (match): $prop_match_min\n";
	print "D\tmaximal percentage of the query size (match): $prop_match_max\n";
	print "Select the parameter you want to change or type Y\n";
	$response=<STDIN>;
	chomp($response);
}
print "running...\n";


$pseudo=0;

foreach(@candidate)
{	
	@linecandidate=split(/\t/, $_);
	chomp(@linecandidate);
	
	if($linecandidate[3] > $linecandidate[2])
	{
		$start_m= $linecandidate[2];
		$end_m = $linecandidate[3];
		$signe = 0;
		$strand = "+";
	}
	else
	{
		$start_m=$linecandidate[3];
		$end_m =$linecandidate[2];
		$signe ++;
		$strand = "-";
	}

	$lgth_m = $end_m +1 - $start_m ;
	$lgth_query=$linecandidate[5];

	$test=0;
	$pseudo=0;
	$stop=0;
	foreach my $feat_object ($seq_object->get_SeqFeatures)
	{
		if($feat_object->primary_tag eq "CDS")
		{
			#retrieve position of the CDS
			$start=$feat_object->start;
			$end=$feat_object->end;
			$strand_g1=$feat_object->strand;
			if($strand_g1 == 1)
			{
				$strand_g = "+";	
			}
			else
			{
				$strand_g = "-";
			}
		}
		
		elsif($feat_object->primary_tag eq "gene")
		{
		#print "dans gene\n";
			#retrieve position of an annotated pseudogene
			if($feat_object->has_tag ("pseudo"))
			{
				$start=$feat_object->start;
				$end=$feat_object->end;
				$strand_g1=$feat_object->strand;
				if($strand_g1 == 1)
				{
					$strand_g = "+";	
				}
				else
				{
					$strand_g = "-";
				}
				
			}
			elsif($feat_object->has_tag ("note"))
			{
				foreach my $dd ($feat_object->get_tag_values("note"))
				{
					if($dd =~ /pseudogene/)
					{
						$start=$feat_object->start;
						$end=$feat_object->end;
						$strand_g1=$feat_object->strand;
						if($strand_g1 == 1)
						{
							$strand_g = "+";	
						}
						else
						{
							$strand_g = "-";
						}
					}	
				}
			}
		}
		
		#does the match correspond to an annotated region?
		
		if(($start_m>= $start && $start_m < $end && $end_m > $start && $end_m <= $end) && $strand_g eq $strand|| ($start_m < $start && $end_m > $start && $end_m <= $end && ($end_m - $start) > 0.3*($end - $start)) && $strand_g eq $strand || ($start_m>= $start && $start_m < $end && $end_m > $end && ($end - $start_m) > 0.3*($end - $start)) && $strand_g eq $strand|| ($start-$start_m > 0 && $start-$start_m <= 200 && $end_m-$end>0 && $end_m-$end<=200) && $strand_g eq $strand)
		{#test if bornes of the match are inside an annotated region (pseudo)
			$test++;
			
			$lgth = $end +1 - $start;
			
			#verification if it is an already annotated pseudogene
			if($feat_object->primary_tag eq "gene")
			{
				if($feat_object->has_tag ("pseudo"))
				{
				 # print "tag pseudo\n";
					if($feat_object->has_tag ("note") && $feat_object->has_tag ("gene"))
					{
						printf OUTFILE "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t\tidentified pseudogene\n", $linecandidate[0],$linecandidate[4],$strand,$start_m,$end_m,$start,$end,$feat_object->get_tag_values("gene"),$feat_object->get_tag_values("note");
						$pseudo=1;
					}
					elsif($feat_object->has_tag ("gene"))
					{
						printf OUTFILE "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t\t\tidentified pseudogene\n", $linecandidate[0],$linecandidate[4],$strand,$start_m,$end_m,$start,$end,$feat_object->get_tag_values("gene");
						$pseudo=1;
					}
					else
					{
						printf OUTFILE "%s\t%s\t%s\t%d\t%d\t%d\t%d\t\t\t\tidentified pseudogene\n", $linecandidate[0],$linecandidate[4],$strand,$start_m,$end_m,$start,$end;
						$pseudo=1;
					}
				}

				elsif($feat_object->has_tag ("note"))
				{
				  #print "tag note\n";
					foreach my $dd ($feat_object->get_tag_values("note"))
					{
						if($dd =~ /pseudogene/)
						{
							printf OUTFILE "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t\tidentified pseudogene\n", $linecandidate[0],$linecandidate[4],$strand,$start_m,$end_m,$start,$end,$feat_object->get_tag_values("gene"),$feat_object->get_tag_values("note");
							$pseudo=1;
						}
						#else #test
						#{
						 #   printf OUTFILE "%s\t%s\t\t\t\t\t\t\t\t\tNOT pseudogene\n", $linecandidate[0],$linecandidate[4];
							
						# } 
					}
				}
			}		
			elsif(($feat_object->primary_tag eq "CDS") && $pseudo==0)
			{
			        #print "tag CDS\n";
				if($feat_object->has_tag ("product"))
				{
					foreach my $dd ($feat_object->get_tag_values("product"))
					{
						if($dd =~ /pseudogene/ && $feat_object->has_tag ("gene"))
						{
							printf OUTFILE "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t\tidentified pseudogene\n", $linecandidate[0],$linecandidate[4],$strand,$start_m,$end_m,$start,$end,$feat_object->get_tag_values("gene"),$feat_object->get_tag_values("product");
							$pseudo=1;
						}
						#else #test
						#{
						#	printf OUTFILE "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t\tNOT pseudogene\n", $linecandidate[0],$linecandidate[4],$strand,$start_m,$end_m,$start,$end,"\t",$feat_object->get_tag_values("product");
						#	$pseudo=1;
  
						#}
					}
				}
			}
		}#test if bornes of the match are inside an annotated region (pseudo)
		
		if((($start_m>= $start && $start_m < $end && $end_m > $start && $end_m <= $end) && $strand_g eq $strand || ($start_m < $start && $end_m > $start && $end_m <= $end && ($end_m - $start) > 0.3*($end - $start)) && $strand_g eq $strand || ($start_m>= $start && $start_m < $end && $end_m > $end && ($end - $start_m) > 0.3*($end - $start)) && $strand_g eq $strand||($start-$start_m > 0 && $start-$start_m <= 200 && $end_m-$end>0 && $end_m-$end<=200)) && $pseudo==0 && $strand_g eq $strand)
		{#test if bornes of the match are inside an annotated region
			$test++;
			$lgth = $end +1 - $start;
			
			#verification of the length of the annotated region compared to the query sequence	
			#if about same length ; test if lgth falls in a distribution restrainted around lgth_query
			if($lgth <= ($lgth_query + 3) && $lgth >= ($lgth_query - 3) && $linecandidate[4]<=$evalue)
			{				
				#verification if it is an already annotated pseudogene
				#verification of the presence of a stop codon
				$cmd="blastdbcmd -db $bank -out temp -entry \"$linecandidate[1]\" -range $start_m-$end_m";
				#printf "$cmd\n";
				system($cmd);
				if($signe != 0) 
				{#if sequence on minus strand		
					$fileb = "temp" ;
 					open(FICHIERB, "<$fileb");
	    				@seq=<FICHIERB>;
					close(FICHIERB);
			
 					$new="";
					foreach(@seq)
					{
						if($_ =~ />/)
						{
							$nom_seq=$_;
						}
						else
						{
							chomp($_);
							$new .= $_;
						}
					}
					$new=reverse($new);
					$new=~ tr/ATCG/TAGC/;
					$taille_seq=length($new);

					#search for a codon stop in the sequence
					for($i=0; $i<($taille_seq-3);)
					{
						$codon=substr($new,$i,3);

						if($codon eq "TAA" || $codon eq "TAG" || $codon eq "TGA")
						{
							$stop++;
							if($feat_object->primary_tag eq "CDS")
							{

								printf OUTFILE "%s\t%s\t%s\t%d\t%d\t%d\t%d\t",$linecandidate[0],$linecandidate[4],$strand,$start_m,$end_m,$start,$end;
								$ratio=$lgth/$lgth_query;
								if($feat_object->has_tag ("gene"))
								{
	 								printf OUTFILE "%s\t", $feat_object->get_tag_values("gene");
								}
								else
								{
									printf OUTFILE "\t";
								}
								if($feat_object->has_tag ("product"))
								{
	 								printf OUTFILE "%s\t", $feat_object->get_tag_values("product");
								}
								else
								{
									printf OUTFILE "\t";
								}					
								if($feat_object->has_tag ("protein_id"))
								{
	 								printf OUTFILE "%s\t", $feat_object->get_tag_values("protein_id");
								}
								else
								{
									printf OUTFILE "\t";
								}					
								printf OUTFILE "bad annotation: presence internal stop codon\t$ratio\n";
							}
							$i=$taille_seq;
						}
						else
						{
							$i=$i+3;
						}	
					}
				}#if sequence on minus strand
				else
				{#if sequence on major strand				
					$fileb = "temp" ;
 					open(FICHIERB, "<$fileb");
					@seq=<FICHIERB>;
					close(FICHIERB);
					$strand="+";
			
					$new="";
					foreach(@seq)
					{	
						
						if($_ =~ />/)
						{
							$nom_seq=$_;
						}
						else
						{
							chomp($_);
							$new .= $_;
						}
					}
			
					#search for a codon stop in the sequence
					$taille_seq=length($new);
			
					for($i=0; $i<($taille_seq-3);)
					{
						$codon=substr($new,$i,3);

						if($codon eq "TAA" || $codon eq "TAG" || $codon eq "TGA")
						{
							$stop++;
							if($feat_object->primary_tag eq "CDS")
							{

								printf OUTFILE "%s\t%s\t%s\t%d\t%d\t%d\t%d\t",$linecandidate[0],$linecandidate[4],$strand,$start_m,$end_m,$start,$end;
								$ratio=$lgth/$lgth_query;
								if($feat_object->has_tag ("gene"))
								{
	 								printf OUTFILE "%s\t", $feat_object->get_tag_values("gene");
								}
								else
								{
									printf OUTFILE "\t";
								}
								if($feat_object->has_tag ("product"))
								{
	 								printf OUTFILE "%s\t", $feat_object->get_tag_values("product");
								}
								else
								{
									printf OUTFILE "\t";
								}					
								if($feat_object->has_tag ("protein_id"))
								{
	 								printf OUTFILE "%s\t", $feat_object->get_tag_values("protein_id");
								}
								else
								{
									printf OUTFILE "\t";
								}					
								printf OUTFILE "bad annotation: presence internal stop codon\t$ratio\n";
							}
							$i=$taille_seq;
						}
						else
						{
							$i=$i+3;
						}	
					}
				}#if sequence on major strand
			}
			
			#else if length annotated region << length query sequence
				#bad annotation of the target genome 			
			elsif($lgth <= ($prop_target*$lgth_query) && $linecandidate[4]<=$evalue)
			{#test length query
				$ratio=$lgth/$lgth_query;
				$cmd="blastdbcmd -db $bank -out temp -entry \"$linecandidate[1]\" -range $start_m-$end_m";
				#printf "$cmd\n";
				system($cmd);
				if($signe != 0) 
				{#if sequence on minus strand		
					$fileb = "temp" ;
 					open(FICHIERB, "<$fileb");
	    				@seq=<FICHIERB>;
					close(FICHIERB);
			
 					$new="";
					foreach(@seq)
					{
						if($_ =~ />/)
						{
							$nom_seq=$_;
						}
						else
						{
							chomp($_);
							$new .= $_;
						}
					}
					$new=reverse($new);
					$new=~ tr/ATCG/TAGC/;
					$taille_seq=length($new);

					#search for a codon stop in the sequence

					for($i=0; $i<($taille_seq-3);)
					{
						$codon=substr($new,$i,3);

						if($codon eq "TAA" || $codon eq "TAG" || $codon eq "TGA")
						{
							$stop++;
							if($feat_object->primary_tag eq "CDS")
							{
								printf OUTFILE "%s\t%s\t%s\t%d\t%d\t%d\t%d\t",$linecandidate[0],$linecandidate[4],$strand,$start_m,$end_m,$start,$end;
								$ratio=$lgth/$lgth_query;
								if($feat_object->has_tag ("gene"))
								{
	 								printf OUTFILE "%s\t", $feat_object->get_tag_values("gene");
								}
								else
								{
									printf OUTFILE "\t";
								}
								if($feat_object->has_tag ("product"))
								{
	 								printf OUTFILE "%s\t", $feat_object->get_tag_values("product");
								}
								else
								{
									printf OUTFILE "\t";
								}					
								if($feat_object->has_tag ("protein_id"))
								{
	 								printf OUTFILE "%s\t", $feat_object->get_tag_values("protein_id");
								}
								else
								{
									printf OUTFILE "\t";
								}					
								printf OUTFILE "bad annotation: short sequence in target genome; internal stop codon\t$ratio\n";
							}
							$i=$taille_seq;
						}
						else
						{
							$i=$i+3;
						}						
					}
				}#if sequence on minus strand
		
				else 
				{# match on the major strand				
					$fileb = "temp" ;
 					open(FICHIERB, "<$fileb");
					@seq=<FICHIERB>;
					close(FICHIERB);
					$strand="+";
			
					$new="";
					foreach(@seq)
					{	
						
						if($_ =~ />/)
						{
							$nom_seq=$_;
						}
						else
						{
							chomp($_);
							$new .= $_;
						}
					}
			
					#search for a codon stop in the sequence
					$taille_seq=length($new);
			
					for($i=0; $i<($taille_seq-3);)
					{
						$codon=substr($new,$i,3);

						if($codon eq "TAA" || $codon eq "TAG" || $codon eq "TGA")
						{
							$stop++;
							if($feat_object->primary_tag eq "CDS")
							{

								printf OUTFILE "%s\t%s\t%s\t%d\t%d\t%d\t%d\t",$linecandidate[0],$linecandidate[4],$strand,$start_m,$end_m,$start,$end;
								$ratio=$lgth/$lgth_query;
								if($feat_object->has_tag ("gene"))
								{
	 								printf OUTFILE "%s\t", $feat_object->get_tag_values("gene");
								}
								else
								{
									printf OUTFILE "\t";
								}
								if($feat_object->has_tag ("product"))
								{
	 								printf OUTFILE "%s\t", $feat_object->get_tag_values("product");
								}
								else
								{
									printf OUTFILE "\t";
								}					
								if($feat_object->has_tag ("protein_id"))
								{
	 								printf OUTFILE "%s\t", $feat_object->get_tag_values("protein_id");
								}
								else
								{
									printf OUTFILE "\t";
								}					
								printf OUTFILE "bad annotation: short sequence in target genome; internal stop codon\t$ratio\n";
							}
							$i=$taille_seq;
						}
						else
						{
							$i=$i+3;
						}	
					}
				}# match on the major strand
					
				if($stop==0) # no stop codon in the match
				{#test case where no stop codon
					if($lgth_m < ($prop_match_min* $lgth_query))
					{#test length match 30%
						if($feat_object->primary_tag eq "CDS")
						{

							printf OUTFILE "%s\t%s\t%s\t%d\t%d\t%d\t%d\t",$linecandidate[0],$linecandidate[4],$strand,$start_m,$end_m,$start,$end;
							$ratio=$lgth/$lgth_query;
							if($feat_object->has_tag ("gene"))
							{
	 							printf OUTFILE "%s\t", $feat_object->get_tag_values("gene");
							}
							else
							{
								printf OUTFILE "\t";
							}
							if($feat_object->has_tag ("product"))
							{
	 							printf OUTFILE "%s\t", $feat_object->get_tag_values("product");
							}
							else
							{
								printf OUTFILE "\t";
							}					
							if($feat_object->has_tag ("protein_id"))
							{
	 							printf OUTFILE "%s\t", $feat_object->get_tag_values("protein_id");
							}
							else
							{
								printf OUTFILE "\t";
							}					
							printf OUTFILE "in CDS but very short match\t$ratio\n";					
						}
					}#test length match 30%
					elsif($lgth_m <=($prop_match_max*$lgth_query))
					{#test length match 80%
						if($feat_object->primary_tag eq "CDS")
						{
							printf OUTFILE "%s\t%s\t%s\t%d\t%d\t%d\t%d\t",$linecandidate[0],$linecandidate[4],$strand,$start_m,$end_m,$start,$end;
							$ratio=$lgth/$lgth_query;
							if($feat_object->has_tag ("gene"))
							{
	 							printf OUTFILE "%s\t", $feat_object->get_tag_values("gene");
							}
							else
							{
									printf OUTFILE "\t";
							}
							if($feat_object->has_tag ("product"))
							{
	 							printf OUTFILE "%s\t", $feat_object->get_tag_values("product");
							}
							else
							{
								printf OUTFILE "\t";
							}					
							if($feat_object->has_tag ("protein_id"))
							{
	 							printf OUTFILE "%s\t", $feat_object->get_tag_values("protein_id");
							}
							else
							{
								printf OUTFILE "\t";
							}					
							printf OUTFILE "bad annotation: short sequence in target genome and in match\t$ratio\n";					
						}
					}#test length match 80%
				}#test case where no stop codon
			}#test length query
			
			#else if length annotated region >> length query sequence
				#probable bad annotation in the query genome
			#elsif($lgth > $lgth_query && $linecandidate[4]<=$evalue)
			#{#test length target
				#$ratio=$lgth/$lgth_query;
#
			#	if($feat_object->primary_tag eq "CDS")
			#	{
			#		printf OUTFILE "%s\t%s\t%s\t%d\t%d\t%d\t%d\t",$linecandidate[0],$linecandidate[4],$strand,$start_m,$end_m,$start,$end;
#
			#		if($feat_object->has_tag ("gene"))
			#		{
	 		#			printf OUTFILE "%s\t", $feat_object->get_tag_values("gene");
			#		}
			#		else
			#		{
			#			printf OUTFILE "\t";
			#		}
			#		if($feat_object->has_tag ("product"))
			#		{
	 		#			printf OUTFILE "%s\t", $feat_object->get_tag_values("product");
			#		}
			#		else
			#		{
			#			printf OUTFILE "\t";
			#		}					
			#		if($feat_object->has_tag ("protein_id"))
			#		{
	 		#			printf OUTFILE "%s\t", $feat_object->get_tag_values("protein_id");
			#		}
			#		else
			#		{
			#			printf OUTFILE "\t";
			#		}					
			#		printf OUTFILE "short sequence in query genome; perhaps pseudogene in query\t$ratio\n";
			#	}
			#}#test length target
			elsif($lgth >($prop_match_max*$lgth_query) && $lgth <$lgth_query && $linecandidate[4]<=$evalue)
			{#length target > 80% length query
				$ratio=$lgth/$lgth_query;
				$cmd="blastdbcmd -db $bank -out temp -entry \"$linecandidate[1]\" -range $start_m-$end_m";
				#printf "$cmd\n";
				system($cmd);
				
				if($signe != 0) #if sequence on minus strand
				{		
					$fileb = "temp" ;
 					open(FICHIERB, "<$fileb");
	    				@seq=<FICHIERB>;
					close(FICHIERB);
			
 					$new="";
					foreach(@seq)
					{
						if($_ =~ />/)
						{
							$nom_seq=$_;
						}
						else
						{
							chomp($_);
							$new .= $_;
						}
					}
					$new=reverse($new);
					$new=~ tr/ATCG/TAGC/;
					$taille_seq=length($new);

					#search for a codon stop in the sequence

					for($i=0; $i<($taille_seq-3);)
					{
						$codon=substr($new,$i,3);

						if($codon eq "TAA" || $codon eq "TAG" || $codon eq "TGA")
						{
							$stop++;
							if($feat_object->primary_tag eq "CDS")
							{
								printf OUTFILE "%s\t%s\t%s\t%d\t%d\t%d\t%d\t",$linecandidate[0],$linecandidate[4],$strand,$start_m,$end_m,$start,$end;
								
								if($feat_object->has_tag ("gene"))
								{
	 								printf OUTFILE "%s\t", $feat_object->get_tag_values("gene");
								}
								else
								{
									printf OUTFILE "\t";
								}
								if($feat_object->has_tag ("product"))
								{
	 								printf OUTFILE "%s\t", $feat_object->get_tag_values("product");
								}
								else
								{
									printf OUTFILE "\t";
								}					
								if($feat_object->has_tag ("protein_id"))
								{
	 								printf OUTFILE "%s\t", $feat_object->get_tag_values("protein_id");
								}
								else
								{
									printf OUTFILE "\t";
								}					
								printf OUTFILE "bad annotation: short sequence in target genome (but >$prop_match_max x query); internal stop codon\t$ratio\n";
							}
							$i=$taille_seq;
						}
						else
						{
							$i=$i+3;
						}	
					}
				}	
				else #match on the major strand
				{				
					$fileb = "temp" ;
 					open(FICHIERB, "<$fileb");
					@seq=<FICHIERB>;
					close(FICHIERB);
					$strand="+";
			
					$new="";
					foreach(@seq)
					{	
						
						if($_ =~ />/)
						{
							$nom_seq=$_;
						}
						else
						{
							chomp($_);
							$new .= $_;
						}
					}
			
					#search for a codon stop in the sequence
					$taille_seq=length($new);
			
					for($i=0; $i<($taille_seq-3);)
					{
						$codon=substr($new,$i,3);

						if($codon eq "TAA" || $codon eq "TAG" || $codon eq "TGA")
						{
							$stop++;
							if($feat_object->primary_tag eq "CDS")
							{

								printf OUTFILE "%s\t%s\t%s\t%d\t%d\t%d\t%d\t",$linecandidate[0],$linecandidate[4],$strand,$start_m,$end_m,$start,$end;

								if($feat_object->has_tag ("gene"))
								{
	 								printf OUTFILE "%s\t", $feat_object->get_tag_values("gene");
								}
								else
								{
									printf OUTFILE "\t";
								}
								if($feat_object->has_tag ("product"))
								{
	 								printf OUTFILE "%s\t", $feat_object->get_tag_values("product");
								}
								else
								{
									printf OUTFILE "\t";
								}					
								if($feat_object->has_tag ("protein_id"))
								{
	 								printf OUTFILE "%s\t", $feat_object->get_tag_values("protein_id");
								}
								else
								{
									printf OUTFILE "\t";
								}					
								printf OUTFILE "bad annotation: short sequence in target genome (but >$prop_match_max x query); internal stop codon\t$ratio\n";
							}
							$i=$taille_seq;
						}
						else
						{
							$i=$i+3;
						}	
					}
				}	
				if($stop == 0)
				{#no stop codon in the match	
					if($feat_object->primary_tag eq "CDS")
					{
						printf OUTFILE "%s\t%s\t%s\t%d\t%d\t%d\t%d\t",$linecandidate[0],$linecandidate[4],$strand,$start_m,$end_m,$start,$end;
						if($feat_object->has_tag ("gene"))
						{
	 						printf OUTFILE "%s\t", $feat_object->get_tag_values("gene");
						}
						else
						{
							printf OUTFILE "\t";
						}					
						if($feat_object->has_tag ("product"))
						{
	 						printf OUTFILE "%s\t", $feat_object->get_tag_values("product");
						}
						else
						{
							printf OUTFILE "\t";
						}					
						if($feat_object->has_tag ("protein_id"))
						{
	 						printf OUTFILE "%s\t", $feat_object->get_tag_values("protein_id");
						}
						else
						{
							printf OUTFILE "\t";
						}					
					
						printf OUTFILE "short sequence in target genome but probably not pseudogene\t$ratio\n";
					}
				}#no stop codon in the match						
			}#length target > 80% length query
		}#test if bornes of the match are inside an annotated region		
	}
	#if nothing matches in the target genome: intergenic
	if ($test==0)
	{
		if($lgth_m > ($prop_match_min* $lgth_query))
		{#test length match
			$cmd="blastdbcmd -db $bank -out temp -entry \"$linecandidate[1]\" -range $start_m-$end_m";
			#printf "$cmd\n";
			system($cmd);
			if($signe != 0) 
			{#if sequence on minus strand		
				$fileb = "temp" ;
 				open(FICHIERB, "<$fileb");
	    			@seq=<FICHIERB>;
				close(FICHIERB);
			
 				$new="";
				foreach(@seq)
				{
			
					if($_ =~ />/)
					{
						$nom_seq=$_;
					}
					else
					{
						chomp($_);
						$new .= $_;
					}
				}
				$new=reverse($new);
				$new=~ tr/ATCG/TAGC/;
				$taille_seq=length($new);
				#print "$taille_seq\n";
				#search for a codon stop in the sequence

				for($i=0; $i<($taille_seq-3);)
				{
					$codon=substr($new,$i,3);

					if(($codon eq "TAA" || $codon eq "TAG" || $codon eq "TGA") && $linecandidate[4]<=$evalue)
					{
						$stop++;
						$i=$taille_seq;
						printf OUTFILE "$linecandidate[0]\t$linecandidate[4]\t$strand\t$start_m\t$end_m\t\t\t\t\t\tintergenic with internal stop codon\n";
					}
					else
					{
						$i=$i+3;
					}	
	
				}
			}#if sequence on minus strand
			else
			{#if sequence on major strand				
				$fileb = "temp" ;
 				open(FICHIERB, "<$fileb");
				@seq=<FICHIERB>;
				close(FICHIERB);
				$strand="+";
			
				$new="";
				foreach(@seq)
				{	
						
					if($_ =~ />/)
					{
						$nom_seq=$_;
					}
					else
					{
						chomp($_);
						$new .= $_;
					}
				}
			
				#search for a codon stop in the sequence
				$taille_seq=length($new);
				#print "$taille_seq\n";
				for($i=0; $i<($taille_seq-3);)
				{
					$codon=substr($new,$i,3);

					if(($codon eq "TAA" || $codon eq "TAG" || $codon eq "TGA") && $linecandidate[4]<=$evalue)
					{
						$stop++;
						printf OUTFILE "$linecandidate[0]\t$linecandidate[4]\t$strand\t$start_m\t$end_m\t\t\t\t\t\tintergenic with internal stop codon\n";	    
						$i=$taille_seq;
					}
					else
					{
						$i=$i+3;
					}	
				}
			}#if sequence on major strand	
			if($stop==0)
			{#no stop codon
				if($lgth_m<=($prop_match_max*$lgth_query) && $linecandidate[4]<=$evalue)
				{
					printf OUTFILE "$linecandidate[0]\t$linecandidate[4]\t$strand\t$start_m\t$end_m\t\t\t\t\t\tintergenic\n";			
				}
				elsif($linecandidate[4]<=$evalue)
				{
					printf OUTFILE "$linecandidate[0]\t$linecandidate[4]\t$strand\t$start_m\t$end_m\t\t\t\t\t\tintergenic but probably not pseudogene\n";			
				}	
			}#no stop codon
		}#test length match
		elsif($linecandidate[4]<=$evalue)
		{
			printf OUTFILE "$linecandidate[0]\t$linecandidate[4]\t$strand\t$start_m\t$end_m\t\t\t\t\t\tintergenic but short match ($lgth_m / $lgth_query)\n";					
		}	
	}
}

close(OUTFILE);
exit;
