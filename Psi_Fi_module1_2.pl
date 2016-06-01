#!/usr/bin/perl
#program which retrieves particular data from a blast file
#Psi-Fi module 1 - Lerat and Ochman (2004) Genome Res.
#version 1_2; last modified 24 june 2014

use Bio::SeqIO;

print "************************ Psi_Fi_module1_2 (version 0.2) **********************\n";

print "Enter the name of the tblastn file (format m=8)\n";

$file = <STDIN>;
chomp($file);
open (FICHIER, "<$file") || die "file not found";
@blast=<FICHIER>;
$nb_ligne=$#blast+1;
print "$nb_ligne\n\n";

print "Enter the name of the target species (as it appears in the tblastn file)\n";
$species=<STDIN> ;
chomp($species);

print "\nEnter the name of the query genome file (genbank format)\n";
$fileb = <STDIN>;
chomp($fileb);

my $seqio_object = Bio::SeqIO->new(-format => 'genbank', -file => "$fileb");
my $seq_object = $seqio_object->next_seq;

print "\n!!! WARNING: the protein ids in the tblastn file (first column) must be the same as the protein ids in the genbank file (/protein_id=) !!!\n\n";

print "Enter the name of the output file\n";
$outfile = <STDIN>;
chomp($outfile);
open (OUTFILE, ">$outfile");

$id_min=79;
$id_max=100;
$dist_merge=300;

print "********************* Parameters *******************\n";
print "A\tminimal percentage of identity: $id_min\n";
print "B\tmaximal percentage of identity: $id_max\n";
print "C\tdistance to merge: $dist_merge\n\n";
print "Select the parameter you want to change or type Y\n\n";
$response=<STDIN>;
chomp($response);

while($response ne "Y")
{
  if($response eq "A")
  {
    print "Enter new value for minimal percentage of identity\n";
    $id_min=<STDIN>;
    chomp($id_min);
  }
  elsif($response eq "B")
  {
    print "Enter new value for maximal percentage of identity\n";
    $id_max=<STDIN>;
    chomp($id_max);
  }
  elsif($response eq "C")
  {
    print "Enter new value for merging\n";
    $dist_merge=<STDIN>;
    chomp($dist_merge);
  }
  else
  {
    print "character not valid\n";
  }
  print "********************* Parameters *******************\n";
  print "A\tminimal percentage of identity: $id_min\n";
  print "B\tmaximal percentage of identity: $id_max\n";
  print "C\tdistance to merge: $dist_merge\n\n";
  print "Select the parameter you want to change or type Y\n\n";
  $response=<STDIN>;
  chomp($response);
}
print "running...\n\n";

$previous="bb";
$n=0;
$ligne=0;
$l=0;

foreach(@blast)
{
  @lineblast=split(/\t/, $_);
  chomp($lineblast[$#lineblast]);
  $l++;

  if($lineblast[1] =~ /$species/)
  {
    if($lineblast[2]>$id_min && $lineblast[2]<=$id_max)
    {
      if($previous ne $lineblast[0])
      {
        $previous = $lineblast[0];
        #$nb_query ++;

        if($ligne != 0)
        {
          printf OUTFILE "$nom[$n]\t$species[$n]\t$position[$n][0]\t$position[$n][1]\t$position[$n][2]\t$taille[$n]\n";
        }

        if($n!=0)
        {
          sub by_begin
          {
            $aa=$$a[0];
            $bb=$$b[0];
            return $aa <=> $bb;
          }
          @tri_position=sort by_begin @position;

          $j=0;
          $m=0;

          for($i=0; $i<$n; $i++)
          {
            $dir=$tri_position[$i][1] - $tri_position[$i][0] ;
            if($dir > 0)
            {
              $diff=abs($tri_position[$i][1]-$tri_position[$i+1][0]);
              if($diff<=$dist_merge)
              {
    #the matches are merged
    $position_finale[$j][0]=$tri_position[$i][0];
    $position_finale[$j][1]=$tri_position[$i+1][1];
    if($tri_position[$i+1][2] <= $tri_position[$i][2])
    {
      $position_finale[$j][2]=$tri_position[$i+1][2];
    }
    else
    {
      $position_finale[$j][2]=$tri_position[$i][2];
    }
    $m++;
              }
              else
              {
    if($m !=0)
    {
      $position_finale[$j+1][0]=$tri_position[$i+1][0];
      $position_finale[$j+1][1]=$tri_position[$i+1][1];
      $position_finale[$j+1][2]=$tri_position[$i+1][2];
      $j++;
      $m=0;
    }
    else
    {
      $position_finale[$j][0]=$tri_position[$i][0];
      $position_finale[$j][1]=$tri_position[$i][1];
      $position_finale[$j][2]=$tri_position[$i][2];
      $j ++;
      if($i==$n-1)
      {
        $position_finale[$j][0]=$tri_position[$i+1][0];
        $position_finale[$j][1]=$tri_position[$i+1][1];
        $position_finale[$j][2]=$tri_position[$i+1][2];
      }
    }
              }
            }
            else
            {
              $diff=abs($tri_position[$i+1][1]-$tri_position[$i][0]);
              if($diff<=$dist_merge)
              {
    $position_finale[$j][0]=$tri_position[$i+1][0];
    $position_finale[$j][1]=$tri_position[$i][1];

    if($tri_position[$i+1][2] <= $tri_position[$i][2])
    {
      $position_finale[$j][2]=$tri_position[$i+1][2];
    }
    else
    {
      $position_finale[$j][2]=$tri_position[$i][2];
    }
    $m++;
              }

              else
              {
    if($m !=0)
    {
      $position_finale[$j+1][0]=$tri_position[$i+1][0];
      $position_finale[$j+1][1]=$tri_position[$i+1][1];
      $position_finale[$j+1][2]=$tri_position[$i+1][2];
      $m=0;
      $j++;
    }
    else
    {
      $position_finale[$j][0]=$tri_position[$i][0];
      $position_finale[$j][1]=$tri_position[$i][1];
      $position_finale[$j][2]=$tri_position[$i][2];
      $j++;
      if($i == $n-1)
      {
        $position_finale[$j][0]=$tri_position[$i+1][0];
        $position_finale[$j][1]=$tri_position[$i+1][1];
        $position_finale[$j][2]=$tri_position[$i+1][2];
      }
    }
              }
            }
          }

          for($i=0; $i<$j+1; $i ++)
          {
            printf OUTFILE "$nom[$i]\t$species[$n]\t$position_finale[$i][0]\t$position_finale[$i][1]\t$position_finale[$i][2]\t$taille[$n]\n";
          }

          @position=();
        }

        $n=0;
        $position[$n][0]=$lineblast[8];
        $position[$n][1]=$lineblast[9];
        $position[$n][2]=$lineblast[10];
        $nom[$n]=$lineblast[0];
        $species[$n]=$lineblast[1];
        $protid=0;

        foreach my $feat_object ($seq_object->get_SeqFeatures)
        {
          if($feat_object->primary_tag eq "CDS")
          {
            if($feat_object->has_tag ("protein_id"))
            {
              foreach my $machin ($feat_object->get_tag_values("protein_id"))
              {
    if($lineblast[0] =~ /$machin/)
    {
      $start=$feat_object->start;
      $end=$feat_object->end;
      $taille[$n]=$end+1-$start;
      $protid=1;
    }
    elsif($protid==0)
    {
      $taille[$n]="NA";
    }
              }
            }
          }
        }
        $ligne++;
      }

      else
      {
        $n ++;
        $position[$n][0]=$lineblast[8];
        $position[$n][1]=$lineblast[9];
        $position[$n][2]=$lineblast[10];
        $nom[$n]=$lineblast[0];
        $species[$n]=$lineblast[1];
        $protid=0;

        foreach my $feat_object ($seq_object->get_SeqFeatures)
        {
          if($feat_object->primary_tag eq "CDS")
          {
            if($feat_object->has_tag ("protein_id"))
            {
              foreach my $machin ($feat_object->get_tag_values("protein_id"))
              {
    if($lineblast[0] =~ /$machin/)
    {
      $start=$feat_object->start;
      $end=$feat_object->end;
      $taille[$n]=$end+1-$start;
      $protid=1;
    }
    elsif($protid==0)
    {
      $taille[$n]="NA";
    }
              }
            }
          }
        }
        $ligne=0;
      }
    }
  }

  if($l == $nb_ligne)
  {

    if($n != 0)
    {
      sub by_begin
      {
        $aa=$$a[0];
        $bb=$$b[0];
        return $aa <=> $bb;
      }
      @tri_position=sort by_begin @position;

      $j=0;
      $m=0;

      for($i=0; $i<$n; $i++)
      {
        $dir=$tri_position[$i][1] - $tri_position[$i][0] ;
        if($dir > 0)
        {
          $diff=abs($tri_position[$i][1]-$tri_position[$i+1][0]);
          if($diff<=$dist_merge)
          {
            $position_finale[$j][0]=$tri_position[$i][0];
            $position_finale[$j][1]=$tri_position[$i+1][1];

            if($tri_position[$i+1][2] <= $tri_position[$i][2])
            {
              $position_finale[$j][2]=$tri_position[$i+1][2];
            }
            else
            {
              $position_finale[$j][2]=$tri_position[$i][2];
            }
            $m++;
          }

          else
          {
            if($m !=0)
            {
              $position_finale[$j+1][0]=$tri_position[$i+1][0];
              $position_finale[$j+1][1]=$tri_position[$i+1][1];
              $position_finale[$j+1][2]=$tri_position[$i+1][2];
              $j++;
              $m=0;
            }
            else
            {
              $position_finale[$j][0]=$tri_position[$i][0];
              $position_finale[$j][1]=$tri_position[$i][1];
              $position_finale[$j][2]=$tri_position[$i][2];
              $j ++;
              if($i==$n-1)
              {
    $position_finale[$j][0]=$tri_position[$i+1][0];
    $position_finale[$j][1]=$tri_position[$i+1][1];
    $position_finale[$j][2]=$tri_position[$i+1][2];
              }
            }
          }
        }
        else
        {
          $diff=abs($tri_position[$i+1][1]-$tri_position[$i][0]);
          if($diff<=$dist_merge)
          {
            $position_finale[$j][0]=$tri_position[$i+1][0];
            $position_finale[$j][1]=$tri_position[$i][1];
            if($tri_position[$i+1][2] <= $tri_position[$i][2])
            {
              $position_finale[$j][2]=$tri_position[$i+1][2];
            }
            else
            {
              $position_finale[$j][2]=$tri_position[$i][2];
            }
            $m++;
          }
          else
          {
            if($m !=0)
            {
              $position_finale[$j+1][0]=$tri_position[$i+1][0];
              $position_finale[$j+1][1]=$tri_position[$i+1][1];
              $position_finale[$j+1][2]=$tri_position[$i+1][2];
              $m=0;
              $j++;
            }
            else
            {
              $position_finale[$j][0]=$tri_position[$i][0];
              $position_finale[$j][1]=$tri_position[$i][1];
              $position_finale[$j][2]=$tri_position[$i][2];
              $j++;
              if($i == $n-1)
              {
    $position_finale[$j][0]=$tri_position[$i+1][0];
    $position_finale[$j][1]=$tri_position[$i+1][1];
    $position_finale[$j][2]=$tri_position[$i+1][2];
              }
            }
          }
        }
      }

      for($i=0; $i<$j+1; $i ++)
      {
        printf OUTFILE "$nom[$i]\t$species[$n]\t$position_finale[$i][0]\t$position_finale[$i][1]\t$position_finale[$i][2]\t$taille[$n]\n";
      }

      @position=();
    }
    else
    {
        printf OUTFILE "$nom[$n]\t$species[$n]\t$position[$n][0]\t$position[$n][1]\t$position[$n][2]\t$taille[$n]\n";		
    }
  }
}
close(FICHIER);
close(OUTFILE);
exit;
