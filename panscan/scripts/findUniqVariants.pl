#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use YAML::XS 'LoadFile';
use perlModules::panscan qw(cleanUp print_help2);
use Cwd;
my $currpath = getcwd();
use lib '$currpath';

# Variable Declaration #
my $inputVcf = 'NA';
my $varType = 'NA';
my $db = 'ALL';
my $overlap = 80;
my $config_file = 'NA';
my $help = undef;
my %types = ('SNP','SNP','INDEL','INDEL','SV','SV');

GetOptions(
	'i=s'   => \$inputVcf,
	't=s'   => \$varType,
	'db=s'  => \$db,
	'op=i'  => \$overlap,
	'config=s'  => \$config_file,
	'help'      => \$help    
) or die "Error in command line arguments\n";

# Specify the path to the config file

if($config_file eq 'NA')
{
	$config_file = 'config.yaml';
}

# Create a new Config::IniFiles object
my $cfg = LoadFile($config_file);

# Check if the config file was successfully read
if(!$cfg)
{
    die("Failed to read YAML file: $config_file\n");
}


# Print the help message #
if ($help)
{
    print_help2();
    exit;
}

$varType = uc($varType);
if($varType eq 'NA')
{
	print "\nERROR : Please provide variant type to compare(SNP/INDEL/SV) !! \n\n";
	print_help2();
	exit;
}
elsif(!defined($types{$varType}))
{
	print "\nERROR : Unknown variant type <$varType>.Please provide variant type to compare(SNP/INDEL/SV) !! \n\n";
	print_help2();
	exit;
}


# Result folders #
my $tmpfolder = 'tmp';
my $resultfolder = "$varType\_Comparison\_Result";

# generating resultfolders #
unless(-d $tmpfolder)
{
	mkdir $tmpfolder;
}

if(-d $resultfolder)
{
	&cleanUp($resultfolder);
	mkdir $resultfolder;
}
else
{
	mkdir $resultfolder;
}

if($varType eq 'SNP')
{
	# Calling subroutene to generate input files for SNP #
	my $inputFile = &generateInput($inputVcf,$tmpfolder,$varType);

	my $sample = 'Input';
	my %DBs = ('dbSNP',$cfg->{databases}->{dbSNP},'gnomAD',$cfg->{databases}->{gnomAD},'1000Genomes',$cfg->{databases}->{Genomes_1000},'GME',$cfg->{databases}->{GME});
	my @db_order = ('dbSNP','gnomAD','1000Genomes','GME');
	my $totalDB = 0;
	if($db ne 'ALL')
	{
		@db_order = ();
		my %selectedDB = ();
		foreach my $eachInpDB(split',',$db)
		{
			if(defined($DBs{$eachInpDB}))
			{
				push(@db_order,$eachInpDB);
				$selectedDB{$eachInpDB} = $DBs{$eachInpDB};
			}
			else
			{
				print "Error : No database available for <$eachInpDB>\n";
				print "Please see the available $varType databases below : \n\n";
				print_help2();
				exit;
			}
		}
		%DBs = %selectedDB;
		$totalDB = scalar(keys %DBs);
	}
	
	my $dbNo = 1;
	foreach my $eachDB(@db_order)
	{
		&compareSNPDB(\$inputFile,$DBs{$eachDB},$sample,$eachDB,$resultfolder,$dbNo,$totalDB);
		$dbNo++;
	}
	print "The SNP comparison completed successfully and the results were generated in the folder : $resultfolder\n";
}
elsif($varType eq 'INDEL')
{
	# Calling subroutene to generate input files for INDEL #
	my $inputFile = &generateInput($inputVcf,$tmpfolder,$varType);

	my $sample = 'Input';
	my %DBs = ('gnomAD',$cfg->{databases}->{GNOMAD_INDEL},'1000Genomes',$cfg->{databases}->{Genomes_1000_INDEL},'GME',$cfg->{databases}->{GME_INDEL});
	my @db_order = ('gnomAD','1000Genomes','GME');
	my $totalDB = 0;
	if($db ne 'ALL')
	{
		@db_order = ();
		my %selectedDB = ();
		foreach my $eachInpDB(split',',$db)
		{
			if(defined($DBs{$eachInpDB}))
			{
				push(@db_order,$eachInpDB);
				$selectedDB{$eachInpDB} = $DBs{$eachInpDB};
			}
			else
			{
				print "Error : No database available for <$eachInpDB>\n";
				print "Please see the available $varType databases below : \n\n";
				print_help2();
				&cleanUp($tmpfolder);
				&cleanUp($resultfolder);
				exit;
			}
		}
		%DBs = %selectedDB;
		$totalDB = scalar(keys %DBs);
	}
	
	my $dbNo = 1;
	foreach my $eachDB(@db_order)
	{
		&compareINDELDB(\$inputFile,$DBs{$eachDB},$sample,$eachDB,$resultfolder,$dbNo,$totalDB);
		$dbNo++;
	}
	print "The INDEL comparison completed successfully and the results were generated in the folder : $resultfolder\n";
}
elsif($varType eq 'SV')
{
	 my %DBs = ('DGV',$cfg->{databases}->{DGV_SV_DEL},'1000Genomes',$cfg->{databases}->{Genomes_1000_SV_DEL});
	my @db_order = ('DGV','1000Genomes');
	my %selectedDB = ();
	if($db ne 'ALL')
	{
		@db_order = ();
                foreach my $eachInpDB(split',',$db)
                {
			if(defined($DBs{$eachInpDB}))
			{
				push(@db_order,$eachInpDB);
				$selectedDB{$eachInpDB} = $DBs{$eachInpDB};
			}
			else
			{
                                print "Error : No database available for <$eachInpDB>\n";
                                print "Please see the available $varType databases below : \n\n";
                                print_help2();
				&cleanUp($tmpfolder);
				&cleanUp($resultfolder);
                                exit;
                        }
                }
		%DBs = %selectedDB;
	}

	# Calling subroutene to generate input files for SV #
	my ($svCount,$inputINS,$inputDEL,$inputDUP,$inputINV) = &generateInputSV($inputVcf,$tmpfolder,$resultfolder);
	if($svCount==0)
	{
		print "WARNING : There is no SVs present in the given vcf file. Pleas make sure the input is a pbsv or recommended SV vcf file\n";
	}
	else
	{
		if(-e $inputINS)
		{
			&compareSVINS($inputINS,$tmpfolder,$resultfolder,$overlap);
		}

		if(-e $inputDEL)
		{
			&compareSVDEL($inputDEL,$tmpfolder,$resultfolder,$overlap,\@db_order,\%DBs);
		}
	}
}

# Cleaning temporary folder #
&cleanUp($tmpfolder);
exit;


#################### SUBROUTINES ###################
sub generateInput
{
	# Accepting arguments #
	my $infile = shift;
	my $resfolder = shift;
	my $db = shift;

	# Resultfiles #
	my $resultfile = "$resfolder\/$db\.txt";

	# Opening the vcf file for reading #
	if($infile=~/\.gz$/)
	{
		open(IN,"gunzip -c $infile |") or die "can't open $infile for reading";
	}
	else
	{
		open IN,"<$infile" or die "Can't open $infile for reading\n";
	}
	
	open OUT,">$resultfile" or die "Can't open $resultfile for writing\n";
	while(<IN>)
	{
		chomp;
		next if(/^\s*$/); # Skipping blank lines present, if any #
		next if(/^\s*\#+.*/); # Skipping headers #

		my @data = split"\t",$_,-1; # split & extract data #
		my $chr = $data[0];
		my $pos = $data[1];
		my $ref = $data[3];
		my $l1 = length($ref);
		my $alt = $data[4];
		
		foreach my $eachAlt(split'\,',$alt,-1) # Multi allelic variants #
		{
			my $l2 = length($eachAlt);
			my $diff = abs($l2-$l1);

			if( ($l1==1) && ($l2==$l1) ) # SNP #
			{
				if($db eq 'SNP')
				{
					print OUT "$chr\t$pos\t$ref\t$eachAlt\n";
				}
			}
			elsif($diff<50) # INDEL #
			{
				if($db eq 'INDEL')
				{
					print OUT "$chr\t$pos\t$ref\t$eachAlt\n";
				}
			}
			else # SV #
			{
				if($db eq 'SV')
				{
					print OUT "$chr\t$pos\t$ref\t$eachAlt\n";
				}
			}
		}
	}
	close IN; # Closing the filehandler #
	close OUT;
	return $resultfile;
}

sub compareSNPDB
{
	my $refInfile = shift;
	my $db = shift;
	my $sample = shift;
	my $dbName = shift;
	my $resfolder = shift;
	my $dbNo = shift;
	my $totalDB = shift;

	my $resultfile1 = '';
	if($dbNo==$totalDB)
	{
		$resultfile1 = "$resfolder\/$sample\-SNP\-$dbName\_UNIQUE_SNPS.txt";
	}
	else
	{
		$resultfile1 = "$resfolder\/$sample\-SNP\-$dbName.txt";
	}
	my $resultfile2 = "$resfolder\/$sample\_SNP\_Stat.txt";

	my %DB_Variants = ();
	my $total = 0;
	open IN1,"<$$refInfile" or die "Can't open $$refInfile for reading\n";
	while(<IN1>)
	{
		chomp;
		next if(/^\s*$/);
		$DB_Variants{$_} = 0;
		$total++;
	}
	close IN1;

	open IN2,"<$db" or die "Can't open $db for reading\n";
	while(<IN2>)
	{
		chomp;
		next if(/^\s*$/);
		if(defined($DB_Variants{$_}))
		{
			$DB_Variants{$_} = 1;
		}
	}
	close IN2;

	my $common = 0;
	open OUT1,">>$resultfile1" or die "Can't open $resultfile1 for writing\n";
	while(my($k,$v) = each %DB_Variants)
	{
		if($v==0)
		{
			print OUT1 "$k\n";
		}
		else
		{
			$common++;
		}
	}
	close OUT1;

	my $unique = $total-$common;
	open OUT2,">>$resultfile2" or die "Can't open $resultfile2 for writing\n";
	if($dbNo==1)
	{
		print OUT2 "Database\tTotal input SNPs\tCommon SNPs\tUnique SNPs\n";
	}
	print OUT2 "$dbName\t$total\t$common\t$unique\n";
	close OUT2;
	
	$$refInfile = $resultfile1;
}

sub compareINDELDB
{
	my $refInfile = shift;
	my $db = shift;
	my $sample = shift;
	my $dbName = shift;
	my $resfolder = shift;
	my $dbNo = shift;
	my $totalDB = shift;

        my $resultfile1 = '';
        if($dbNo==$totalDB)
        {
                $resultfile1 = "$resfolder\/$sample\-SNP\-$dbName\_UNIQUE_INDELS.txt";
        }
        else
        {
		$resultfile1 = "$resfolder\/$sample\_INDEL\-$dbName.txt";
        }

	my $resultfile2 = "$resfolder\/$sample\_INDEL\_Stat.txt";

	my %DB_Variants = ();
	my $total = 0;
	open IN1,"<$$refInfile" or die "Can't open $$refInfile for reading\n";
	while(<IN1>)
	{
		chomp;
		next if(/^\s*$/);
		$DB_Variants{$_} = 0;
		$total++;
	}
	close IN1;

	open IN2,"<$db" or die "Can't open $db for reading\n";
	while(<IN2>)
	{
		chomp;
		next if(/^\s*$/);
		if(defined($DB_Variants{$_}))
		{
			$DB_Variants{$_} = 1;
		}
	}
	close IN2;

	my $common = 0;
	open OUT1,">>$resultfile1" or die "Can't open $resultfile1 for writing\n";
	while(my($k,$v) = each %DB_Variants)
	{
		if($v==0)
		{
			print OUT1 "$k\n";
		}
		else
		{
			$common++;
		}
	}
	close OUT1;

	my $unique = $total-$common;
	open OUT2,">>$resultfile2" or die "Can't open $resultfile2 for writing\n";
	if($dbNo==1)
	{
		print OUT2 "Database\tTotal input INDELs\tCommon INDELs\tUnique INDELs\n";
	}
	print OUT2 "$dbName\t$total\t$common\t$unique\n";
	close OUT2;
	
	$$refInfile = $resultfile1;
}

sub generateInputSV
{
	my $infile = shift;
	my $resFolder = shift;
	my $resFolder2 = shift;
	my $sample = 'Input';

	my $resultfile1 = "$resFolder\/$sample\_SV_INS.txt";
	my $resultfile2 = "$resFolder\/$sample\_SV_DEL.txt";
	my $resultfile3 = "$resFolder\/$sample\_SV_DUP.txt";
	my $resultfile4 = "$resFolder\/$sample\_SV_INV.txt";
	my $resultfile5 = "$resFolder2\/$sample\_SV_Stat.txt";

	my %VarInfo = ();
	my @header = ();
	open OUT1,">$resultfile1" or die "Can't open $resultfile1 for writing\n";
	open OUT2,">$resultfile2" or die "Can't open $resultfile2 for writing\n";
	open OUT3,">$resultfile3" or die "Can't open $resultfile3 for writing\n";
	open OUT4,">$resultfile4" or die "Can't open $resultfile4 for writing\n";
	open IN,"<$infile" or die "Can't open $infile for reading\n";
	while(<IN>)
	{
		chomp;
		next if(/^\s*$/);
		if(/^\s*\#CHR.*/)
		{
			@header = split'\t',$_,-1;
		}
		next if(/^\s*\#+.*/);
		my @data = split"\t",$_,-1;
		if($data[6] eq 'PASS')
		{
			my ($SVType,$SVLen,$end) = ('NA','NA','NA');
			foreach my $eachInfo(split'\;',$data[7],-1)
			{
				my($k,$v) = split'\=',$eachInfo,-1;
				if($k eq 'SVTYPE')
				{
					$SVType = $v;
				}
				elsif($k eq 'END')
				{
					$end = $v;
				}
				elsif($k eq 'SVLEN')
				{
					$SVLen = $v;
				}
			}
			my $chr = $data[0];
			$chr=~s/chr//gi;
			if($chr eq 'X')
			{
				$chr = 23;
			}
			if($chr eq 'Y')
			{
				$chr = 24;
			}
			
			if($SVType eq 'INS')
			{
				$end = $data[1]+$SVLen;
				if(length($chr)==1)
				{
					print OUT1 "$chr\t$data[1]\t$end\n";
				}
			}
			elsif($SVType eq 'DEL')
			{
				if(length($chr)==1)
				{
					print OUT2 "$chr\t$data[1]\t$end\n";
				}
			}
			elsif($SVType eq 'DUP')
			{
				print OUT3 "$chr\t$data[1]\t$end\n";
			}
			elsif($SVType eq 'INV')
			{
				print OUT4 "$chr\t$data[1]\t$end\n";
			}
			
			foreach my $i(9..$#header)
			{
				my $gt = (split'\:',$data[$i],-1)[0];
				if( ($gt eq '1/0') ||($gt eq '1/1') ||($gt eq '0/1') ||($gt eq './1') || ($gt eq '1/.'))
				{
					$VarInfo{$header[$i]}{'TYPE'}{$SVType}++;
					$VarInfo{$header[$i]}{'TYPE'}{'TOTAL'}++;
				}
			}
		}
	}
	close IN;
	close OUT1;
	close OUT2;
	close OUT3;
	close OUT4;

	open OUT5,">$resultfile5" or die "Can't open $resultfile5 for writing\n";
	my @SVTypes = ('INS','DEL','DUP','INV','BND');
	print OUT5 "Sample\tTotal Variants\t",(join"\t",@SVTypes)."\n";
	my $totSVs = 0;
	foreach my $i(9..$#header)
	{
		my $sample = $header[$i];
		my @sample_result = ($sample);
		foreach my $eachType('TOTAL',@SVTypes)
		{
			if(defined($VarInfo{$header[$i]}{'TYPE'}{$eachType}))
			{
				$totSVs++;
				push(@sample_result,$VarInfo{$header[$i]}{'TYPE'}{$eachType});
			}
			else
			{
				push(@sample_result,0);
			}
		}
		my $eachResult = join"\t",@sample_result;
		print OUT5 "$eachResult\n";
	}
	close OUT5;
	if(-z $resultfile1)
	{
		unlink $resultfile1;
	}
	if(-z $resultfile2)
	{
		unlink $resultfile2;
	}
	if(-z $resultfile3)
	{
		unlink $resultfile3;
	}
	if(-z $resultfile4)
	{
		unlink $resultfile4;
	}
	return ($totSVs,$resultfile1,$resultfile2,$resultfile3,$resultfile4);
}

sub compareSVINS
{
	my $infile = shift;
	my $tmpfolder = shift;
	my $resultfolder = shift;
	my $svOverlap = shift;

	my $resultfile = "$tmpfolder\/After-comparison-with-DGV_Final-Uniq-SVs.txt";
	my $retval = system("java $cfg->{classes}->{reciprocal_overlap} $infile $cfg->{databases}->{DGV_SV_INS} $svOverlap > $resultfile");
	if($retval==0)
	{
		&parseSVINS($resultfile,$resultfolder);
		print "SV insertion comparison completed !!\n";
	}
	else
	{
		print "Problem in the reciprocal overlap\n";
		exit;
	}
}

sub parseSVINS
{
	my $infile = shift;
	my $resultfolder = shift;
	my $percentage = shift;

	my $resultfile1 = "$resultfolder\/After-comparison-with-DGV_Final-Uniq-SV-Insertions.txt";
	my $resultfile2 = "$resultfolder\/SV-Insetion_Comparison_Stat.txt";

	my $total = 0;
	my $common = 0;
	my $unique = 0;

	open OUT1,">$resultfile1" or die "Can't open $resultfile1 for writng\n";
	open IN,"<$infile" or die "Can't open $infile for reading\n";
	while(<IN>)
	{
		chomp;
		next if(/^\s*$/);
		$total++;
		my @data = split"\t",$_,-1;

		if($data[1]>0)
		{
			$common++;
		}
		elsif($data[1]==0)
		{
			my($chr,$cord) = split'\:',$data[0],-1;
			my($cord1,$cord2) = split'\-',$cord,-1;
			$chr=~tr/Chr//d;
			print OUT1 "$chr\t$cord1\t$cord2\n";
			$unique++;
		}
	}
	close IN;
	close OUT1;

	open OUT2,">$resultfile2" or die "Can't open $resultfile2 for writng\n";
	print OUT2 "Total insertions\tCommon insertions\tUnique insertions\n";
	print OUT2 "$total\t$common\t$unique\n";
	close OUT2;
}

sub compareSVDEL
{
	my $infile = shift;
	my $tmpfolder = shift;
	my $resultfolder = shift;
	my $svOverlap = shift;
	my $refDBorder = shift;
	my $refDBs = shift;

	my $totalDBs = scalar(@$refDBorder);
	my $no = 1;
	my $retFile = '';
	foreach my $eachDb(@$refDBorder)
	{
		my $dbfile = $$refDBs{$eachDb};
		my $resultfile = '';
		$resultfile = "$tmpfolder\/After-comparison-with\-$eachDb\-SV-Deletions.txt";
		if($no==$totalDBs)
		{
			$resultfile = "$tmpfolder\/After-comparison-with\-$eachDb\_Final-Uniq-SV-Deletions.txt";
		}
		
		my $retval = 'NA';

		if($no==1)
		{
			$retval = system("java $cfg->{classes}->{reciprocal_overlap} $infile $dbfile $svOverlap > $resultfile");
		}
		else
		{
			$retval = system("java $cfg->{classes}->{reciprocal_overlap} $retFile $dbfile $svOverlap > $resultfile");
		}
		if($retval==0)
		{
			$retFile = &parseSVDEL($resultfile,$resultfolder,$eachDb,$totalDBs,$no);
			if($no==$totalDBs)
			{
				print "SV Deletion comparison completed !!\n";
			}
		}
		else
		{
			print "Problem in the reciprocal overlap\n";
			exit;
		}
		$no++;
	}
}

sub parseSVDEL
{
	my $infile = shift;
	my $resultfolder = shift;
	my $db = shift;
	my $totalDB = shift;
	my $N = shift;

	my $resultfile1 = "$resultfolder\/After-comparison-with\-$db\_SV-Deletions.txt";
	if($N==$totalDB)
	{
		$resultfile1 = "$resultfolder\/After-comparison-with\-$db\_Final-Uniq-SV-Deletions.txt";
	}
	my $resultfile2 = "$resultfolder\/SV-Deletion_Comparison_Stat.txt";

	my $total = 0;
	my $common = 0;
	my $unique = 0;

	open OUT1,">$resultfile1" or die "Can't open $resultfile1 for writing\n";
	open IN,"<$infile" or die "Can't open $infile for reading\n";
	while(<IN>)
	{
		chomp;
		next if(/^\s*$/);
		$total++;
		my @data = split"\t",$_,-1;

		if($data[1]>0)
		{
			$common++;
		}
		elsif($data[1]==0)
		{
			my($chr,$cord) = split'\:',$data[0],-1;
			my($cord1,$cord2) = split'\-',$cord,-1;
			$chr=~tr/Chr//d;
			print OUT1 "$chr\t$cord1\t$cord2\n";
			$unique++;
		}
	}
	close IN;
	close OUT1;

	if($N==1)
	{
		open OUT2,">$resultfile2" or die "Can't open $resultfile2 for writng\n";
		print OUT2 "Database\tTotal deletions\tCommon deletions\tUnique deletions\n";
	}
	else
	{
		open OUT2,">>$resultfile2" or die "Can't open $resultfile2 for writng\n";
	}

	print OUT2 "$db\t$total\t$common\t$unique\n";
	close OUT2;

	return $resultfile1;
}
