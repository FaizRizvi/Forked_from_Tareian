#!/usr/bin/perl 
use strict;
no warnings 'uninitialized';


if ((@ARGV==0) || ($ARGV[0] eq "-h"))
 	{
 		usage();
 		die"\n";
 	}
	
if ($ARGV[0] eq "-v")
	{
		die "\nPeaKDEck v1.1; Feb 2014\n\n";
	}


my $genomeSizeFile = 0;
my $genomeBaseCount = 0;
my @lineArray = ();
my $lineCount = 0;
my $totalArea = 0;
my $count = 0;
my $currentStep = 0;
my @wordArray = ();

use POSIX 'acos'; 
use POSIX 'fmod';
#use Math::Random::MT qw(srand rand); <-windows only

#stamp package:
use Cwd;
my $dir = getcwd;
my $user = getlogin();
#my @userData = getpwnam($user);
#my $hd = $#userData - 1;
my @timeData = split(" ", localtime(time));
my $timeString = join(" ", @timeData);
#my $shell = $#userData;
my $progLoc = $0;
my $command = join(" ", @ARGV);
my $sys = $^O;


print STDERR "\nProgram: PeaKDEck ($progLoc)\nVersion: v1.1\n\nstart time: $timeString\nuser: $user\ncwd: $dir\nparameters: $command\nsystem: $sys\n\n";
#end stamp package

##default values for density analyzer and peak caller, subject to change by user command line input
my $stepSize = 50;

my $peakBinSize = 300;
my $backgroundBinSize = 3000;
my $tail = $peakBinSize/2;

our $sigma = 50;
my $thresh = "DEFAULT";
my $offSet = 0;
my $maxThresh = 100000000;
my $bandWidth = $tail*2 +1;
my $binSize = $tail*2 +1;

my $thresholdSignificanceLevel = 0.001;
my $peaksForBackgroundCalc = 50000;

our %densityFrequencyHash;
my %pValueTable = ();
my $pcrDuplicates = "OFF";

my $closedAndOpen = "NO";
my $maxPeakPValue = "OFF";
my $mapqThresh = 0;
my $uqThresh = 10000;

my $peakStringency = 1;
my $powerDigit = 5;

my $numberDataSetsPerPeak = 60;

#end default values



#####################
my $peakThresholdValue = 0;
#####################

# check to ensure that an even number of items exist on the command line before converting @ARGV to a hash
die "ERROR: uneven option-value pairings.\n\n" if @ARGV%2 != 0;
my %arguments = @ARGV;
my @arguments = ();

#### checks to make sure that only one function is selected
my $processCount = 0;
$processCount += 1 if exists $arguments{"-D"};
$processCount += 1 if exists $arguments{"-F"};
$processCount += 1 if exists $arguments{"-P"};
$processCount += 1 if exists $arguments{"-R"};

$processCount += 1 if exists $arguments{"-T"};
$processCount += 1 if exists $arguments{"-NS"};


die "too many processes selected.\n\n" if $processCount > 1;
die "no processes selected.\n\n" if $processCount < 1;
############################################################

##### runs the random read sam selection option
my @randomTitleCheck = ();
my $samForRandom = 0;
my $randomNumber = 0;

if (exists $arguments{"-R"})
	{
		#check correct file extension
		
		$samForRandom = $arguments{"-R"}; delete($arguments{"-R"});

		die "Input file (set to \'$samForRandom\') not found.\n\n" unless -e $samForRandom;
		my @samTestResults = fileTest('SAM', $samForRandom);
		print STDERR "\n\tINPUT file format test report \'SAM\': $samForRandom\n\n\t\tSAM format ... $samTestResults[0]\n\t\tSAM read order ... $samTestResults[1]\n\t\tSAM chr grouping ... $samTestResults[3]\n\t\tSAM header ... $samTestResults[2]\n\n";
		die "sam file format for inputfile (\'$samForRandom\') not recognised.\n\n"	if $samTestResults[0] eq 'BAD';
		print STDERR "warning: sam file reads for inputfile (\'$samForRandom\') appear disordered.\n\n"	if $samTestResults[1] eq 'BAD';
		print STDERR "warning: sam file chromosomes for inputfile (\'$samForRandom\') appear disordered.\n\n"	if $samTestResults[3] eq 'BAD';

			
		#assign variables and delete from %arguments hash
		$randomNumber = $arguments{"-nr"} if exists $arguments{"-nr"}; delete($arguments{"-nr"});
				
		#checks for and kills if extra command line arguments
		@arguments = %arguments;
		die "unused command line arguments: @arguments\n\n" if keys %arguments > 0;
		
		#sets limits for assigned parameters
		die "random number of reads ($randomNumber; -nr) must be a positive integer. \n\n" unless (($randomNumber =~ /[0-9]/) && (fmod($randomNumber,1) == 0) && ($randomNumber > 0)); 
		
		#executes subroutine	 
		print STDERR "Random Read Selector v1.0\n\nwith parameters:\n\tinput file: $samForRandom\n\ttarget number of reads: $randomNumber\n\n";
		randomReadSelection($samForRandom, $randomNumber);# <- need to write subroutine
	}

##### runs the top peaks selection option
my @topPeaksTitleCheck = ();
my $bedForTop = 0;
my $topPeaksNumber = 'ALL';

if (exists $arguments{"-T"})
	{
		$bedForTop = $arguments{"-T"}; delete($arguments{"-T"});
		#check correct file extension
		die "Input file (set to \'$bedForTop\') not found.\n\n" unless -e $bedForTop;

		my @bedTestResults = fileTest('BED', $bedForTop);
		print STDERR "\n\tINPUT file format test report \'BED\': $bedForTop\n\n\t\tBED format ... $bedTestResults[0]\n\n";
		die "bed file format for inputfile (\'$bedForTop\') not recognised.\n\n"	if $bedTestResults[0] eq 'BAD';
						
		#assign variables and delete from %arguments hash
		$topPeaksNumber = $arguments{"-n"} if exists $arguments{"-n"}; delete($arguments{"-n"});
				
		#checks for and kills if extra command line arguments
		@arguments = %arguments;
		die "unused command line arguments: @arguments\n\n" if keys %arguments > 0;
		
		#sets limits for assigned parameters()
		die "number of peaks (set to \'$topPeaksNumber\') must be either 'ALL' or a positive integer. \n\n" unless ((($topPeaksNumber =~ /[0-9]/) & (fmod($topPeaksNumber,1) == 0) && ($topPeaksNumber >= 0)) ||  ($topPeaksNumber eq 'ALL'));
				
		#executes subroutine	 
		print STDERR "Top Peak Selector v1.0\n\nwith parameters:\n\tinput file: $bedForTop\n\ttarget number of reads: $topPeaksNumber\n\n";
		topPeaks($bedForTop, $topPeaksNumber);
	}

##### runs the sam numerical sort option
my @samNSortTitleCheck = ();
my $samForNSort = 0;

if (exists $arguments{"-NS"})
	{
		$samForNSort = $arguments{"-NS"}; delete($arguments{"-NS"});

		die "Input file (set to \'$samForNSort\') not found.\n\n" unless -e $samForNSort;
		my @samTestResults = fileTest('SAM', $samForNSort);
		print STDERR "\n\tINPUT file format test report \'SAM\': $samForNSort\n\n\t\tSAM format ... $samTestResults[0]\n\t\tSAM read order ... $samTestResults[1]\n\t\tSAM chr grouping ... $samTestResults[3]\n\t\tSAM header ... $samTestResults[2]\n\n";
		die "sam file format for inputfile (\'$samForNSort\') not recognised.\n\n"	if $samTestResults[0] eq 'BAD';
				
		#checks for and kills if extra command line arguments
		@arguments = %arguments;
		die "unused command line arguments: @arguments\n\n" if keys %arguments > 0;
		
		#executes subroutine	 
		print STDERR "Numerical Sam Sort v1.0\n\nwith parameters:\n\tinput file: $samForNSort\n\n";
		samNSort($samForNSort);
	}
##### runs the sam filter option
my @filterTitleCheck = ();
my @headerTitleCheck = ();
my $samForFilter = 0;
my $samForHeader = 0;
if (exists $arguments{"-F"})
	{
		$samForFilter = $arguments{"-F"}; delete($arguments{"-F"});
		die "Input file (set to \'$samForFilter\') not found.\n\n" unless -e $samForFilter;
		my @samTestResults = fileTest('SAM', $samForFilter);
		print STDERR "\n\tINPUT file format test report \'SAM\': $samForFilter\n\n\t\tSAM format ... $samTestResults[0]\n\t\tSAM read order ... $samTestResults[1]\n\t\tSAM chr grouping ... $samTestResults[3]\n\t\tSAM header ... $samTestResults[2]\n\n";
		die "sam file format for inputfile (\'$samForFilter\') not recognised.\n\n"	if $samTestResults[0] eq 'BAD';
		die "sam file reads for inputfile (\'$samForFilter\') appear disordered.\n\n"	if $samTestResults[1] eq 'BAD';
		print STDERR "warning: sam file chromosomes for inputfile (\'$samForFilter\') appear disordered.\n\n"	if $samTestResults[3] eq 'BAD';

		die "genome size file (-g) must be set.\n\n" unless exists $arguments{"-g"};
		
		# checks that -g has been set (mandatory)
		$genomeSizeFile = $arguments{"-g"} if exists $arguments{"-g"}; delete($arguments{"-g"});
		die "chromosome size file (set to \'$genomeSizeFile\') required but not found.\n\n" unless ((defined $genomeSizeFile) && (-e $genomeSizeFile));	
		my @genomeTestResults = fileTest('CHROM', $genomeSizeFile);
		print STDERR "\n\tchromosome size file format test report \'CHROM\': $genomeSizeFile\n\n\t\tCHROM format ... $genomeTestResults[0]\n\n";
		die "chromosome size file format for inputfile (\'$genomeSizeFile\') not recognised.\n\n"	if $genomeTestResults[0] eq 'BAD';

		
			
		#assign targe file name
		#$samForFilter = $arguments{"-F"}; delete($arguments{"-F"});
		$samForHeader = $samForFilter unless exists $arguments{"-i"};
		$samForHeader = $arguments{"-i"} if exists $arguments{"-i"};
		

		if (exists $arguments{"-i"})
			{
				die "$samForHeader (-i) not found.\n\n" unless -e $samForHeader; 
				my @samTestHeaderResults = fileTest('SAM', $samForHeader);
				print STDERR "\n\tINPUT file format test report \'SAM\': $samForHeader\n\n\t\tSAM format ... $samTestHeaderResults[0]\n\t\tSAM read order ... $samTestHeaderResults[1]\n\t\tSAM chr grouping ... $samTestHeaderResults[3]\n\t\tSAM header ... $samTestHeaderResults[2]\n\n";
				die "sam file header for inputfile (\'$samForHeader\') is not recognised.\n\n"	if $samTestHeaderResults[2] eq 'NONE';
				delete($arguments{"-i"});
			}
		
		  		
		#assign variables and delete from %arguments hash
		$mapqThresh = $arguments{"-q"} if exists $arguments{"-q"}; delete($arguments{"-q"});
		$uqThresh = $arguments{"-u"} if exists $arguments{"-u"}; delete($arguments{"-u"});
		$pcrDuplicates = $arguments{"-PCR"} if exists $arguments{"-PCR"}; delete($arguments{"-PCR"});
		
		#checks for and kills if extra command line arguments
		@arguments = %arguments;
		die "unused command line arguments: @arguments\n\n" if keys %arguments > 0;
		
		#sets limits for assigned parameters
		die "mapq threshold (set to \'$mapqThresh\') must be a non-negative integer.\n\n" unless (($mapqThresh =~ /[0-9]/) && (fmod($mapqThresh,1) == 0) && ($mapqThresh >= 0)); 
		die "uq threshold (set to \'$uqThresh\') must be a non-negative integer.\n\n" unless (($uqThresh =~ /[0-9]/) && (fmod($uqThresh,1) == 0) && ($uqThresh >= 0)); 
		
		#die "$samForHeader (-i) not found.\n\n" if !-e $samForHeader; 
		die "PCR duplicates (set to \'$pcrDuplicates\') can only have the values ON (remove PCR duplicates) or OFF (retain PCR duplicates).\n\n" if (($pcrDuplicates ne "OFF") && ($pcrDuplicates ne "ON")); 
		
		#executes subroutine	 
		print STDERR "Sam Filter v1.0.\n\nwith parameters:\n\tinput file: $samForFilter\n\theader file: $samForHeader\n\tmapq threshold: $mapqThresh\n\tuq threshold: $uqThresh\n\tPCR duplicate removal: $pcrDuplicates\n\n";
		samFilter($samForHeader, $samForFilter, $mapqThresh, $uqThresh);
	}

##### runs the density reader option
my $samForDensity = 0;
my @samTitleCheck = ();

if (exists $arguments{"-D"})
	{
		
		$samForDensity = $arguments{"-D"}; delete($arguments{"-D"});
		die "Input file (set to \'$samForDensity\') not found.\n\n" unless -e $samForDensity;
		
		#check correct file type				
		my @samTestResults = fileTest('SAM', $samForDensity);
		print STDERR "\n\tINPUT file format test report \'SAM\': $samForDensity\n\n\t\tSAM format ... $samTestResults[0]\n\t\tSAM read order ... $samTestResults[1]\n\t\tSAM chr grouping ... $samTestResults[3]\n\t\tSAM header ... $samTestResults[2]\n\n";
		die "sam file format for inputfile (\'$samForDensity\') not recognised.\n\n"	if $samTestResults[0] eq 'BAD';
		die "sam file reads for inputfile (\'$samForDensity\') appear disordered.\n\n"	if $samTestResults[1] eq 'BAD';
		die "sam file chromosomes for inputfile (\'$samForDensity\') appear disordered.\n\n"	if $samTestResults[3] eq 'BAD';

		# checks that -g has been set (mandatory)
		die "genome size file (-g) must be set.\n\n" unless exists $arguments{"-g"};
		$genomeSizeFile = $arguments{"-g"} if exists $arguments{"-g"}; delete($arguments{"-g"});
		die "chromosome size file (set to \'$genomeSizeFile\') required but not found.\n\n" unless ((defined $genomeSizeFile) && (-e $genomeSizeFile));	
		my @genomeTestResults = fileTest('CHROM', $genomeSizeFile);
		print STDERR "\n\tchromosome size file format test report \'CHROM\': $genomeSizeFile\n\n\t\tCHROM format ... $genomeTestResults[0]\n\n";
		die "chromosome size file format for inputfile (\'$genomeSizeFile\') not recognised.\n\n"	if $genomeTestResults[0] eq 'BAD';
		
		#assign variables and delete from %arguments hash
		$thresh = $arguments{"-t"} if exists $arguments{"-t"}; delete($arguments{"-t"});
		$maxThresh = $arguments{"-m"} if exists $arguments{"-m"}; delete($arguments{"-m"});
		$tail = $arguments{"-n"} if exists $arguments{"-n"}; delete($arguments{"-n"});
		$sigma = $arguments{"-d"} if exists $arguments{"-d"}; delete($arguments{"-d"});
		$stepSize = $arguments{"-STEP"} if exists $arguments{"-STEP"}; delete($arguments{"-STEP"});
		$offSet = $arguments{"-o"} if exists $arguments{"-o"}; delete($arguments{"-o"});
		$genomeSizeFile = $arguments{"-g"} if exists $arguments{"-g"}; delete($arguments{"-g"});
		
		#checks for and kills if extra command line arguments
		@arguments = %arguments;
		die "unused command line arguments: @arguments\n\n" if keys %arguments > 0;
		
		die "minimum bin threshold (set to \'$thresh\') be 'DEFAULT' or a non-negative integer.\n\n" unless ((($thresh =~ /[0-9]/) && (fmod($thresh,1) == 0) && ($thresh >= 0)) || ($thresh eq 'DEFAULT')); 
		die "maximum bin threshold (set to \'$maxThresh\') be a positive integer.\n\n" unless (($maxThresh =~ /[0-9]/) && (fmod($maxThresh,1) == 0) && ($maxThresh > 0)); 
		die "tailSize must (set to \'$tail\') be a positive integer.\n\n" unless (($tail =~ /[0-9]/) && (fmod($tail,1) == 0) && ($tail > 0)); 
		die "sigma must (set to \'$sigma\') be a positive integer.\n\n" unless (($sigma =~ /[0-9]/) && (fmod($sigma,1) == 0) && ($sigma > 0)); 
		die "step size must (set to \'$stepSize\') be a positive integer.\n\n" unless (($stepSize =~ /[0-9]/) && (fmod($stepSize,1) == 0) && ($stepSize > 0)); 
		die "offset must (set to \'$offSet\') be an integer.\n\n" unless (($offSet =~ /[0-9]/) && (fmod($offSet,1) == 0));
		
		#executes subroutine	 
		print STDERR "Density Analyzer v1.0.\n\nwith parameters: \n\tinput file: $samForDensity\n\tstep size: $stepSize\n\ttail length: $tail\n\tsigma: $sigma\n\tminimum bin threshold: $thresh\n\tmaximum bin threshold: $maxThresh\n\ttrack offset: $offSet\n\tbin size: $binSize\n\n";
		densityAnalyzer($samForDensity, $stepSize, $tail, $sigma, $thresh, $offSet, $maxThresh, $genomeSizeFile);
	}

##### runs the peakCaller option
my $bluePrintFile = "NONE";
my $samForReport = 0;
my $samRPBScore = 0;
@samTitleCheck = ();
my @blueTitleCheck = ();
our $flatThreshold = "NONE";

if (exists $arguments{"-P"})
	{
		#### NB thresholdSignificanceLevel to build in  
		#check correct file extension
		$samForReport = $arguments{"-P"}; delete($arguments{"-P"});
		die "Input file (set to \'$samForReport\') not found.\n\n" unless -e $samForReport;
		my @samTestResults = fileTest('SAM', $samForReport);
		print STDERR "\n\tINPUT file format test report \'SAM\': $samForReport\n\n\t\tSAM format ... $samTestResults[0]\n\t\tSAM read order ... $samTestResults[1]\n\t\tSAM chr grouping ... $samTestResults[3]\n\t\tSAM header ... $samTestResults[2]\n\n";
		die "sam file format for inputfile (\'$samForReport\') not recognised.\n\n"	if $samTestResults[0] eq 'BAD';
		die "sam file reads for inputfile (\'$samForReport\') appear disordered.\n\n"	if $samTestResults[1] eq 'BAD';
		die "sam file chromosomes for inputfile (\'$samForReport\') appear disordered.\n\n"	if $samTestResults[3] eq 'BAD';
		

		$closedAndOpen = "YES" if exists $arguments{"-b"};
		$bluePrintFile = $arguments{"-b"} if exists $arguments{"-b"}; 

		if (exists $arguments{"-b"})
			{
				die "Blueprint file (set to \'$bluePrintFile\') not found.\n\n" unless -e $bluePrintFile;
				my @bedTestResults = fileTest('BED', $bluePrintFile);
				print STDERR "\n\tBLUEPRINT file format test report \'BED\': $bluePrintFile\n\n\t\tBED format ... $bedTestResults[0]\n\n";
				die "bed file format for inputfile (\'$bluePrintFile\') not recognised.\n\n"	if $bedTestResults[0] eq 'BAD';
				delete($arguments{"-b"});
			}
		
		# checks that -g has been set (mandatory)
		die "genome size file (-g) must be set.\n\n" unless exists $arguments{"-g"};
		$genomeSizeFile = $arguments{"-g"} if exists $arguments{"-g"}; delete($arguments{"-g"});
		die "chromosome size file (set to \'$genomeSizeFile\') required but not found.\n\n" unless ((defined $genomeSizeFile) && (-e $genomeSizeFile));	
		my @genomeTestResults = fileTest('CHROM', $genomeSizeFile);
		print STDERR "\n\tchromosome size file format test report \'CHROM\': $genomeSizeFile\n\n\t\tCHROM format ... $genomeTestResults[0]\n\n";
		die "chromosome size file format for inputfile (\'$genomeSizeFile\') not recognised.\n\n"	if $genomeTestResults[0] eq 'BAD';


		#assign variables and delete from %arguments hash
		$peakBinSize = $arguments{"-bin"} if exists $arguments{"-bin"}; delete($arguments{"-bin"});
		$backgroundBinSize = $arguments{"-back"} if exists $arguments{"-back"}; delete($arguments{"-back"});
		$peaksForBackgroundCalc = $arguments{"-npBack"} if exists $arguments{"-npBack"}; delete($arguments{"-npBack"});
		$thresholdSignificanceLevel = $arguments{"-sig"} if exists $arguments{"-sig"}; delete($arguments{"-sig"});
		$maxPeakPValue = $arguments{"-PVAL"} if exists $arguments{"-PVAL"}; delete($arguments{"-PVAL"});
		$flatThreshold = $arguments{"-FLAT"} if exists $arguments{"-FLAT"}; delete($arguments{"-FLAT"});
		$stepSize = $arguments{"-STEP"} if exists $arguments{"-STEP"}; delete($arguments{"-STEP"});
				
		#checks for and kills if extra command line arguments
		@arguments = %arguments;
		die "unused command line arguments: @arguments\n\n" if keys %arguments > 0;
		
		#sets limits for assigned parameters
		die "peak bin size (set to \'$peakBinSize\') must be a positive integer.\n\n" unless (($peakBinSize !~ /[^0-9\.]/) && (fmod($peakBinSize,1) == 0) && ($peakBinSize > 0)); 
		die "background bin size (set to \'$backgroundBinSize\') must be a positive integer.\n\n" unless (($backgroundBinSize !~ /[^0-9\.]/) && (fmod($backgroundBinSize,1) == 0) && ($backgroundBinSize > 0));  
		die "step size (set to \'$stepSize\') must be a positive integer.\n\n" unless (($stepSize !~ /[^0-9\.]/) && (fmod($stepSize,1) == 0) && ($stepSize > 0)); 
		die "number of peaks for background sampling (set to \'$peaksForBackgroundCalc\') must be a positive integer.\n\n" unless (($peaksForBackgroundCalc !~ /[^0-9\.]/) && (fmod($peaksForBackgroundCalc,1) == 0) && ($peaksForBackgroundCalc > 0)); 
		die "flat threshold (set to \'$flatThreshold\') must be 'NONE' or a non-negative number.\n\n" unless ((($flatThreshold !~ /[^0-9\.]/) && ($flatThreshold >= 0)) || ($flatThreshold eq 'NONE')); 
		die "threshold significance level (set to \'$thresholdSignificanceLevel\') must be in the range [0 .. 1].\n\n" unless (($thresholdSignificanceLevel >= 0) && ($thresholdSignificanceLevel <= 1)); 

		die "background bin size (set to \'$backgroundBinSize\') must be greater than peak bin size (set to \'$peakBinSize\').\n\n" if $backgroundBinSize < $peakBinSize; 
		die "peak score p values (-PVAL; set to \'$maxPeakPValue\') must be either ON (for p values) or OFF (for maximum density score)\n\n" if (($maxPeakPValue ne "OFF") && ($maxPeakPValue ne "ON"));								

		#executes subroutine	 
		print STDERR "Peak Caller 1.1\n\nwith parameters:\n\tinput sam file: $samForReport\n\tblueprint file: $bluePrintFile\n\tbin size/background size: $peakBinSize/$backgroundBinSize\n\tpeak target for threshold/binscore assessment: $peaksForBackgroundCalc\n\tcalculated peak threshold significance: $thresholdSignificanceLevel\n\tstep size: $stepSize\n\tflat threshold: $flatThreshold\n\tpeak score p values: $maxPeakPValue\n\n";
		my @thresholdValue = threshold($samForReport, $bluePrintFile, $peakBinSize, $backgroundBinSize, $peaksForBackgroundCalc, $genomeSizeFile, $thresholdSignificanceLevel, $closedAndOpen, $flatThreshold, $maxPeakPValue); #if $flatThreshold eq "NONE";### <- need to add powerDigit to subroutine
		$thresholdValue[0] = $flatThreshold if $flatThreshold ne 'NONE';
		peakCall($samForReport, $thresholdValue[0], $peakBinSize, $backgroundBinSize, $genomeSizeFile, $stepSize, $thresholdValue[1], $maxPeakPValue);
	}


my $RSname = 0;
my $RSnumber = 0;


###########################################################################
############################## sub routines ###############################
###########################################################################


sub readCount
	{
		my $RCtarget = $_[0];
		my $RCcount = 0;
		my $buffer = 0;
		my @RCarray = ();
		my $RCheaderCount = 0;
		open (TARGET, "<", $RCtarget);
		while (<TARGET>)
			{
				@RCarray = split ("", $_);
				$RCheaderCount += 1 if $RCarray[0] eq "@";
				last if $RCarray[0] ne "@";
			}
		close TARGET;
		
		open (my $targetHandle, "<", $RCtarget);
		while(sysread $targetHandle, $buffer, 4096)
			{
				$RCcount += ($buffer =~ tr/\n//);
			}

		close $targetHandle;
		
		$RCcount = $RCcount - $RCheaderCount;
		return $RCcount;
	}

# sub randomSam
# 	{
# 		my $outputfile = $_[0];
# 		$outputfile = "$outputfile" . ".sam";
# 		my $targetReadNumber = $_[1];
# 		my $chromosomeSizeFile = $_[2];

# 		my $index = 0;
# 		my $x = 0;
# 		my $range = 0;
# 		my @line = ();
# 		my $key = 0;
# 		my $readName = 0;
# 		my $randomReadNumber = 0;
# 		my @chromName = ();
# 		my @chromSize = ();
# 		my $chromCount = 0;
# 		my $end = 0;
# 		my $start = 0;
# 		my $total = 0;
# 		my $readStart = 0;
# 		my @chromStart = ();
# 		my @chromEnd = ();
# 		my $chrom = 0;
# 		my %indexHash = ();
# 		my $multipleCount = 0;
		

# 		open (CHROMS, "<", "$chromosomeSizeFile");
# 		while (<CHROMS>)
# 			{	
# 				@line = split(" ", $_);
# 				$end = $line[1] + $total;
# 				$start =  1 + $total; 
# 				push (@chromName, $line[0]);
# 				push (@chromStart, $start);
# 				push (@chromEnd, $end);
# 				$total = $total + $line[1];
# 				$chromCount += 1;
# 			}
# 		close CHROMS;

# 		$range = $total - 1;

# 		while ($randomReadNumber < $targetReadNumber)
# 			{
# 				$index = int(rand($range)) + 1;
# 				$indexHash{$index} += 1;
# 				$randomReadNumber+=1;
# 			}

# 		print STDERR "\t$targetReadNumber reads selected...\n";
# 		$randomReadNumber = 0;
		
# 		open (OUTPUT, ">", $outputfile);
# 		for ($x=0; $x<$total; $x +=1)
# 			{
				
# 				if (exists $indexHash{$x})
# 					{ 
# 				while ($multipleCount < $indexHash{$x})
# 					{
# 						#$randomReadNumber+= 1;
# 						#$readName = "#" ."$randomReadNumber";
# 						print OUTPUT "#rand\t0\t$chrom\t$readStart\t37\t37M\t*\t0\t0\tRANDOM_Seq\tRANDOM_Qual\tUQ:i:0\tNM:i:0\n" if exists $indexHash{$x};
# 						$multipleCount += 1;
# 						#print STDERR "1,000,000 reads done\n" if $randomReadNumber%1000000 == 0;
# 					}
# 				$multipleCount = 0;
# 			}
# 				print STDERR "$x bases covered \n" if $x%100000000 == 0;
# 			}

# 		# $randomReadNumber = 0;
# 		# open (OUTPUT, ">", $outputfile);

# 		# for $key ( sort {$a<=>$b} keys %indexHash) 
# 		# 	{
# 		# 		for ($x = 0; $x < $chromCount; $x += 1)
# 		# 			{
# 		# 				while ($multipleCount < $indexHash{$key})
# 		# 					{
# 		# 						if (("$key" >= $chromStart[$x]) && ("$key" <= $chromEnd[$x]))
# 		# 							{
# 		# 								$randomReadNumber+= 1;
# 		# 								$readStart = "$key" - $chromStart[$x];
# 		# 								$chrom = $chromName[$x];
# 		# 								$readName = "#" ."$randomReadNumber";
# 		# 								print OUTPUT "$readName\t0\t$chrom\t$readStart\t37\t37M\t*\t0\t0\tRANDOM_Seq\tRANDOM_Qual\tUQ:i:0\tNM:i:0\n";
# 		# 							}
# 		# 						$multipleCount += 1;
# 		# 					}
# 		# 				$multipleCount = 0;
# 		# 			}
# 		# 	}
# 		# close OUTPUT;
# 		close OUTPUT;
# 		print STDERR "\t$outputfile file created.\n\n";
# 	}


sub threshold 
	{
		my $samForAnalysis = $_[0];
		my $bluePrint = $_[1];
		my $halfRegionSize = ($_[2])/2;
		my $halfBackgroundSize = ($_[3])/2;
		my $Size = $_[2];
		my $peakBinSize = $Size;
		my $numberOfPeaks = $_[4]; # the number of sites used to calculate p value: default 10000; probably won't give option to change this.
		my $genomeSizeFile = $_[5];
		my $thresholdSignificanceLevel = $_[6];
		my $closedAndOpen = $_[7];
		my $flatThreshold = $_[8];
		my $maxPeakPValue =  $_[9];

		my @wholeBlue = ();
		my $index = 0;
		my $open = 0;
		
		my $closed = 0;
		my $x = 0;
		my $range = 0;

		my %closedChrom = ();
		my %openChrom = ();
		
		my @closedChromArray = ();
		my @openChromArray = ();
		
		my @line = ();
		my $siteLine = 0;
		my @globalSelection = ();
		my $key = 0;
		my $count = 0;
		my $subRange = 0;
		my $readName = 0;
		my $thresholdFromThreshold = 0;

		my $start = 0;
		my $end = 0;
		my @chromLine = ();
		my @chromName = ();
		my @chromStart = ();
		my @chromEnd = ();

		my $total = 0;
		my $chromCount = 0;
		
		my $randomPeakNumber = 0;
		my $randomPeakStart = 0;
		my $randomPeakEnd = 0;
		my $randomPeakString = 0;
		my $chrom = 0;
		my @randomChromArray = ();

		my $randomRef = 0;
		
		my %randomChromHash = ();
		my $randomStringData = 0;
		my @random = ();
		
		my $randomReadsPerBin = ();
		my $randomPeaksCompleted = ();

		my $openRef = 0;
		
		my $openStringData = 0;
		my @open = ();
		
		my $openReadsPerBin = ();
		my $openPeaksAssessed = ();

		my $closedRef = 0;
		
		my $closedStringData = 0;
		my @closed = ();
		
		my $closedReadsPerBin = ();
		my $closedPeaksCompleted = ();
		my $pValueTableRef;
		my $numberDataSetsPerPeak = 60;

		my $closedCount = 0;
		my $openCount = 0;
		my $costaCounter = 0;
		my $sum = 0;
		my @blueLine = ();
		
		my @openOnes;
		my @closedOnes;
		
		if ($closedAndOpen eq "YES")
			{
				open (BLUE, "<", "$bluePrint") || die  "Unable to open $bluePrint\n\n";
				print STDERR "\treading blueprint file ($bluePrint)...\n\n";
				
				while (<BLUE>)
					{
						push (@wholeBlue, $_);
						@blueLine = split(" ", $_);
						push (@openOnes, $_) if $blueLine[4] >= $numberDataSetsPerPeak;
						push (@closedOnes, $_) if $blueLine[4] == 0;
					}
				close BLUE;
				
				print STDERR "\trandomly selecting ~$numberOfPeaks open/closed sites...\n\n";
				
				@globalSelection = ();
				$range = @closedOnes;
				@globalSelection = (0) x $range;
				
				if ($range <= $numberOfPeaks)
					{
						@closedChromArray = @closedOnes;
					}
				
				if ($range > $numberOfPeaks)
					{
						while ($closed < $numberOfPeaks) 
							{
								$index = int(rand($range));
								next if $globalSelection[$index] > 0;
								$closedChrom{$index} = $closedOnes[$index];
								$closed += 1;
								$globalSelection[$index] = 1;
							}

						$x = 0;
						for $key ( sort {$a<=>$b} keys %closedChrom) 
							{
						           $closedChromArray[$x] = $closedChrom{$key};
						           $x += 1;
							}
					}
				my $peakTally;
				$peakTally = @closedChromArray;	
				
				@globalSelection = (0) x $range;
				$range = @openOnes;
				$globalSelection[$range] = 0;
				for (@globalSelection)
					{
						$globalSelection[$_] = 0;
					}

				if ($range <= $numberOfPeaks)
					{
						@openChromArray = @openOnes;
					}
				
				if ($range > $numberOfPeaks)
					{
						while ($open < $numberOfPeaks) 
							{
								$index = int(rand($range));
								next if $globalSelection[$index] > 0;
								$openChrom{$index} = $openOnes[$index];
								$open += 1;
								$globalSelection[$index] = 1;
							}

						$x = 0;
						for $key ( sort {$a<=>$b} keys %openChrom) 
							{
						           $openChromArray[$x] = $openChrom{$key};
						           $x += 1;
							}
					}
				
				$peakTally = @openChromArray;	
				
				$openRef = \@openChromArray;
				$closedRef = \@closedChromArray;

				# free some memory! 20130425
				%closedChrom = ();
				%openChrom = ();
				@wholeBlue = ();
 
				# calculate the mean read counts in each site category in the target sam file
				
				print STDERR "\tFor open signal:\n";
				
				
				@open = signalIntensity($openRef, 0, $samForAnalysis, $numberOfPeaks);
				$openReadsPerBin = $open[1];
				$openPeaksAssessed = $open[2];

				print STDERR "$openReadsPerBin reads per bin ($peakBinSize bp) at $openPeaksAssessed randomly selected sites common to $numberDataSetsPerPeak of 125 open chromatin data sets\n\n"; 
				print STDERR "\tFor closed signal:\n";
				
				@closed = signalIntensity($closedRef, 1, $samForAnalysis, $numberOfPeaks);
				
				$thresholdFromThreshold = $closed[0];
				$closedReadsPerBin = $closed[1];
				$closedPeaksCompleted = $closed[2];
				$pValueTableRef = $closed[3];

				print STDERR "\t$closedReadsPerBin reads per bin ($peakBinSize bp) at $closedPeaksCompleted randomly selected ENCODE closed chromatin sites.\n\n";
				print STDERR "\tcalculated read density threshold (p < $thresholdSignificanceLevel): $thresholdFromThreshold\n";
				print STDERR "\tusing flat threshold at $flatThreshold for peak calling.\n\n" if $flatThreshold ne "NONE";
				
			}
		
		if ($closedAndOpen eq "NO")
			{

				### HERE
				print STDERR "\tLooking for chromosomes in $samForAnalysis... \n";
				my $startTime = time;
				my %chromsPresent;
				my $headerOff = 0;
				open (my $targetChroms, "<", $samForAnalysis);
				while (<$targetChroms>)
					{
						
						if ($headerOff == 0)
							{
								my @header = split("", $_);
								$headerOff = 1 if $header[0] ne "@";
								next if $headerOff == 0;
							}
						$_ =~ /(chr\w*)/;
						$chromsPresent{$1} = 1 unless exists $chromsPresent{$1};
					}
				close $targetChroms;
				my $endTime = time;
				my $duration = $endTime - $startTime;
				### TO HERE

				open (CHROMS, "<", "$genomeSizeFile");
				while (<CHROMS>)
					{	
						@chromLine = split(" ", $_);
						$end = $chromLine[1] + $total;
						$start =  1 + $total; 
						push (@chromName, $chromLine[0]) if $chromsPresent{$chromLine[0]} == 1; # HERE
						
						push (@chromStart, $start)if $chromsPresent{$chromLine[0]} == 1; # HERE
						push (@chromEnd, $end)if $chromsPresent{$chromLine[0]} == 1; # HERE
						$total = $total + $chromLine[1];
						$chromCount += 1;
						
					}
				
				die "\tnone of the chromosomes listed in $genomeSizeFile have been detected in $samForAnalysis.\n\n" unless scalar(@chromName) >= 1;
				close CHROMS;

				$range = $total - 1;
				my $randomPeakNumber = 0;
				
				while ($randomPeakNumber < $numberOfPeaks)
					{
						$index = int(rand($range)) + 1;
												
						for ($x = 0; $x < $chromCount; $x += 1)
							{
								if (($index >= $chromStart[$x]) && ($index <= $chromEnd[$x]))
									{
										$randomPeakNumber+= 1;
										$randomPeakStart = $index - $chromStart[$x];
										$randomPeakEnd = $randomPeakStart + 10;
										$chrom = $chromName[$x];
										$readName = "#" ."$randomPeakNumber";
										$randomPeakString = "$chrom\t$randomPeakStart\t$randomPeakEnd\t$readName\n";
										$randomChromHash{$index} = $randomPeakString;
									}	
							} 	
					}	

				#need to sort the randomchromarray 
				$x = 0;
				for $key ( sort {$a<=>$b} keys %randomChromHash) 
					{
				           $randomChromArray[$x] = $randomChromHash{$key};
				           $x += 1;
					}
				$randomRef = \@randomChromArray;
				my @randomStringData = signalIntensity($randomRef, 1, $samForAnalysis, $numberOfPeaks, $genomeSizeFile, $_[2], $_[3], $thresholdSignificanceLevel, $maxPeakPValue);
				@random = split(" ", $randomStringData);

				$thresholdFromThreshold = $randomStringData[0];
				$randomReadsPerBin = $randomStringData[1];
				$randomPeaksCompleted = $randomStringData[2];
				$pValueTableRef = $randomStringData[3];

				print STDERR "\t$randomReadsPerBin reads per bin ($peakBinSize bp) at $randomPeaksCompleted random sites.\n";
				print STDERR "\tcalculated read density threshold (p < $thresholdSignificanceLevel): $thresholdFromThreshold\n";
				print STDERR "\tusing flat threshold at $flatThreshold for peak calling.\n\n" if $flatThreshold ne "NONE";
			}

		return ($thresholdFromThreshold, $pValueTableRef);
	}




# # # # # # # # # # # # # # # # # 
# Sub to measure signal intensity  within each randomly selected target site
# # # # # # # # # # # # # # # # # 

sub signalIntensity
	{
		print STDERR "\tassessing signal intensity...\n\n";
		my @tedPeaksArray = @{$_[0]};
		my @bedPeaksArray = @tedPeaksArray;

		my $kernalOn = $_[1];

		my $numberOfSubPeaks = $_[3];
		my $targetSamFile = $_[2];
		my $sizeOfBin = $peakBinSize;
		my $sizeOfBackground = $backgroundBinSize;
		my $distanceFromPeakCentre = $sizeOfBin/2;
		my $backgroundFromPeakCentre = $sizeOfBackground/2;
				
		my $reCalc = 1;
		my $currentTargetPeak = 0;		
		my $lowPeak = 0;	
		my $highPeak = 0;		
		my $peakCentre = 0;
		my $peakChrom = 0;
		my $peakName = 0;
		my $lowLimit = 0; 
		my $highLimit = 0;
		
		my $actualPeakCount = 0;
		my $currentValue = 0;
		
		my $headerLineCheckerOff = 0;
		my $s = 0;
		my $flyCount = 0;

		my @globalBin = ();
		my @splitArray = ();
		my @headerCheck = ();
				
		my @splitSamArray = ();
		my @lineChromSizes = ();
		my %chromOrderScore = ();
		my $simpleCount = 1;
		my @scoreArray = ();

		my @probabilityArray = ();
		my $backgroundScore = 0;
		my $backgroundLowLimit = 0;
		my $backgroundHighLimit = 0;
		my $correctedFlyCount = 0;
		my $meanBackground = 0;
		my $totalFlyCount = 0;
		
		my $returnString = 0;
		my $thresholdFromSignalIntensity = 0;
		my $peaksCompleted = 0;
		my $arrayChecker = 0;
		my $readsPerBin = 0;

		my @lowBack = ();
		my @highBack = ();
		my @lowBin = ();
		my @highBin = ();

		my @centrePoint = ();
		my $x = 0;
		my @lineSplitter = ();
		my @chromosomeNamesForIntensity = 0;

		# get the order of chromosomes to exclude skips;
		open (CHROM, "<", "$genomeSizeFile") or die "$genomeSizeFile not found\n\n";
		while (<CHROM>)
			{
				@lineChromSizes = split(" ", $_);
				$chromOrderScore{"$lineChromSizes[0]"} = $simpleCount;
				$simpleCount += 1;
			}
		close CHROM;
		

		$actualPeakCount = @bedPeaksArray;
		
		print STDERR "\t$actualPeakCount sites were selected randomly.\n";
		
		for ($x = 0; $x<$actualPeakCount; $x += 1)
			{
				@lineSplitter = split (" ", $bedPeaksArray[$x]);
				$centrePoint[$x] = int(($lineSplitter[1] + $lineSplitter[2])/2);
				$chromosomeNamesForIntensity[$x] = $lineSplitter[0];
			}

# make the outer bins, and exclude any overlapping bins
my $lowBackCount = 0;
my $highBackCount = 0;
my $checkingCounters = 0;
my $nextEx = 0;
my $loopCount = 0;
my $diffCount = 0;
my $spliced = 0;
my $totalRawCentralBinScore = 0;
my $totalCheckingCounters = 0;
while ($loopCount <= 10)
	{		
		$actualPeakCount = @centrePoint;
		
		for ($x = 0; $x<$actualPeakCount; $x += 1)
					{
						$lowBack[$x] = $centrePoint[$x] - $backgroundFromPeakCentre;
						$highBack[$x] = $centrePoint[$x] + $backgroundFromPeakCentre;
					}
		
		$lowBackCount = @lowBack;
		$highBackCount = @highBack;
		
		for ($x = 0; $x<$actualPeakCount; $x += 1)
					{
						$nextEx = $x + 1;
						
						if (($highBack[$x] > $lowBack[$nextEx]) && ($chromosomeNamesForIntensity[$x] eq $chromosomeNamesForIntensity[$nextEx])) 
							{
								$checkingCounters += 1;
								$spliced = splice (@centrePoint, $nextEx, 1);
								splice (@lowBack, $nextEx, 1);
								splice (@highBack, $nextEx, 1);
								splice (@chromosomeNamesForIntensity, $nextEx, 1);
							} 
					}

		$loopCount += 1;
		$totalCheckingCounters = $totalCheckingCounters + $checkingCounters;
		last if $checkingCounters == 0;
		$checkingCounters = 0;
	}

print STDERR "\t$totalCheckingCounters sites discarded for overlapping.\n\n";
# make the boundary coordinates of inner bins
$actualPeakCount = @centrePoint;
for ($x = 0; $x<$actualPeakCount; $x += 1)
		{
			$lowBin[$x] = $centrePoint[$x] - $distanceFromPeakCentre;
			$highBin[$x] = $centrePoint[$x] + $distanceFromPeakCentre;
		}
		
# find out how many reads in each bin:
my $chr = 0;
my $lowBackground = 0;
my $lowCentre = 0;
my $highCentre = 0;
my $highBackground = 0;
my $targetLine = 0;
my @samArray = ();
my $backGroundClock = 0;
my $binClock = 0;
my @centreScoresArray = ();
my @backgroundScoresArray = ();
my $chrNameCheck = 0;
my $loopChromNameCheck = 0;

# From HERE:
# open (SAM, "<", $targetSamFile);
# while (<SAM>)
# 	{
# 		$targetLine = $_;
# 		@samArray = split(" ", $targetLine);
# 		@headerCheck = split("", $_);
# 		last if $headerCheck[0] ne "@";
# 		$headerLineCheckerOff = 1 if $headerCheck[0] ne "@";
# 	}

## check sites on each chromosome
my %chromSiteCount = ();
for ($x=0; $x<$actualPeakCount; $x+=1)
	{
		$chromSiteCount{"$chromosomeNamesForIntensity[$x]"} += 1;
	}

my $lastChrom = 0;
		$headerLineCheckerOff = 0;
		open (SAM, "<", $targetSamFile);
		while (<SAM>)
			{
				if ($headerLineCheckerOff == 0)
					{
						@headerCheck = split("", $_);
						$headerLineCheckerOff = 1 if $headerCheck[0] ne "@";
						$targetLine = $_ if $headerCheck[0] ne "@";
						next if $headerLineCheckerOff == 0;
						last if $headerLineCheckerOff == 1;
					}
			}		

		$targetLine =~ /\b(chr\w*)\b\t\b(\w*)\b/;
			my $lineChrom = $1;
			my $basePos = $2;

		for ($x=0; $x<=$#chromosomeNamesForIntensity; $x+=1)
			{
				if ($lastChrom ne $chromosomeNamesForIntensity[$x])
					{
						print STDERR "\tsampling signal at $chromSiteCount{$chromosomeNamesForIntensity[$x]} sites on $chromosomeNamesForIntensity[$x]...\n";
						
						if ($chromosomeNamesForIntensity[$x] ne $lineChrom)
							{
								while ($chromosomeNamesForIntensity[$x] ne $lineChrom)
									{
										$targetLine = <SAM>;
										$targetLine =~ /\b(chr\w*)\b\t\b(\w*)\b/;
										$lineChrom = $1;
										$basePos = $2;
										last if $targetLine = eof;
									}
							}
						}
				
				while (($lineChrom eq $chromosomeNamesForIntensity[$x]) && ($basePos < $highBack[$x]))
					{
						$backGroundClock += 1 if (($basePos > $lowBack[$x]) && ($basePos < $highBack[$x]));
						$binClock += 1 if (($basePos > $lowBin[$x]) && ($basePos < $highBin[$x]));
						
						$targetLine = <SAM>;
						$targetLine =~ /\b(chr\w*)\b\t\b(\w*)\b/;
						$lineChrom = $1;
						$basePos = $2;
						last if $targetLine = eof;
						#print STDERR "$basePos\t";
					}

				push(@centreScoresArray, $binClock);
				push(@backgroundScoresArray, $backGroundClock);
				
				$binClock = 0;
				$backGroundClock = 0;
				$lastChrom = $chromosomeNamesForIntensity[$x];
			}
# TO HERE

my @standardArray = ();
my $correctedFlyCountPlus5000 = 0;

# take bin scores and make background corrected density scores
my $maximumBinScore = 0;
my $backgroundCorrectedDensity;
#%densityFrequencyHash is global to allow calculation of pvalueTable.

for ($x=0; $x<$actualPeakCount; $x+=1)
	{
		$meanBackground = ($backgroundScoresArray[$x]/$sizeOfBackground)*$sizeOfBin;
		$meanBackground = sprintf("%.0f", $meanBackground); #rounds the number of reads  per bin to integer;
		$meanBackground = int($meanBackground); #converts to int
		$backgroundCorrectedDensity = $centreScoresArray[$x] - $meanBackground;
		$densityFrequencyHash{$backgroundCorrectedDensity} += 1;

		$maximumBinScore = $backgroundCorrectedDensity if $backgroundCorrectedDensity > $maximumBinScore;

		$totalRawCentralBinScore = $totalRawCentralBinScore + $centreScoresArray[$x];
	}
# calculate some numbers		
$peaksCompleted = @centreScoresArray;

$readsPerBin = $totalRawCentralBinScore/$peaksCompleted;
$readsPerBin = sprintf("%.2f", $readsPerBin);

my $standardArrayCount = @standardArray;
		# $threshold and @dataForPValue $readsPerBin and $peaksCompleted are global variables. $threshold is used and must be set here prior to peak calling (peak calling sub will crash without it). @dataForPValue is a list of bin scores in closed chromatin, and can be called on by the peakPValue sub. $readsPerBin and $peaks completed are collected on the fly by threshold sub.

		$thresholdFromSignalIntensity = kernel(\%densityFrequencyHash, $thresholdSignificanceLevel) if $kernalOn == 1;
		#make a hash for looking up peak P values
		my %pValueTable = ();
		if (($kernalOn == 1) && ($maxPeakPValue eq "ON"))
			{
				print STDERR "\n\tcreating peak density p value table...\n";
				
				for ($x=0; $x<=50000; $x+=1)
					{
						$pValueTable{$x} = peakPValue($x);
						$x=50001 if $pValueTable{$x} == "0.000000000001";
					}
			}

		@probabilityArray = ();
		@bedPeaksArray = ();
		my $pValueTableRef = \%pValueTable;
		print STDERR "\n\t$peaksCompleted sites sampled.\n\n";
		
		$returnString = "$thresholdFromSignalIntensity " . "$readsPerBin " . "$peaksCompleted";
		return  ($thresholdFromSignalIntensity, $readsPerBin, $peaksCompleted, \%pValueTable);
	}

# sub kernel runs only once: it takes an array of bin scores from the closed chromatin sites, and returns a threshold corresponding to the globally set $thresholdSignificanceLevel
	
sub kernel 
	{
		my %hashForKernel = %{$_[0]};
		my $thresholdSignificanceLevel = $_[1];
		
		my $n;
		my $sumTotal;
		my $squareDifference;
		my $runningTotal;
		my $difference;

		#calculate standard deviation of values for bandwidth estimation:
		foreach my $key ( keys %hashForKernel )
			{
				$n += $hashForKernel{$key};
			}
		
		my $bandWidth = 1;
		my $oneOverRootTwoPi = 0.39894228;
		my $e = 2.718281828;
		my $round = 0;
		my $probTotal = 0;
		my $eFactor = 0;
		my $eFigure = 0;
		my $gaussFigure = 0;
		my $rankNumber = 0;
		my $fracFactor = 0;
		my $interTotal = 0;
		my $actualTotal = 0;
		my @rankProbArray = ();
		my $posCount = 0;
		my $diff = 0;
		my $absDiff = 0;
		my $thresholdFromKernel = 0;

		

		for ($rankNumber = 0; $rankNumber <= 100; $rankNumber += 0.1)
			{
				
				foreach my $key ( keys %hashForKernel )
					{
					  for (my $x = 0; $x < $hashForKernel{$key}; $x += 1)
					  	{

					  		$diff = $rankNumber - $key;
							$absDiff = abs($diff);
							$eFactor = (($absDiff*$absDiff)/$bandWidth)/-2;
							
							$eFigure = $e**$eFactor;
							$gaussFigure = $oneOverRootTwoPi*$eFigure;
												
							$interTotal = $gaussFigure;
							$probTotal = $probTotal + $interTotal;
					  	}
					}

				$actualTotal = $probTotal/($n*$bandWidth);
				# stops the looping as soon as P is less than $thresh significance level
				if ($actualTotal < $thresholdSignificanceLevel)
					{
						$thresholdFromKernel = sprintf("%.1f", $rankNumber);
						for ($rankNumber = 0; $rankNumber < 100; $rankNumber += 1)
								{
									my $indexing = $rankNumber/10;
								}	
						@rankProbArray = (); # frees memory;
						return $thresholdFromKernel;
						last;
					}

				push (@rankProbArray, $actualTotal);
				$probTotal = 0;
				$posCount = 0;
			}
	}


# calculates p values for individual peak maxima, using the global data array, @dataForPValue, which is filled by the signalIntensity sub.
sub peakPValue
	{
		my $maximaForP = $_[0];
		
		my $nP = 0;

		foreach my $key ( keys %densityFrequencyHash )
			{
				$nP += $densityFrequencyHash{$key};
			}

		die "array with sample bin scores is empty!\n\n" if ($nP == 0);
		
		my $bandWidth = 1;
		my $oneOverRootTwoPi = 0.39894228;
		my $e = 2.718281828;
		my $round = 0;
		my $probTotal = 0;
		my $eFactor = 0;
		my $eFigure = 0;
		my $gaussFigure = 0;
		my $fracFactor = 0;
		my $interTotal = 0;
		my $actualTotal = 0;
		my $diff = 0;
		my $absDiff = 0;
		
		foreach my $key ( keys %densityFrequencyHash )
			{
			  for (my $x = 0; $x < $densityFrequencyHash{$key}; $x += 1)
				  	{

						$diff = $maximaForP - $key;
						$absDiff = abs($diff);
						$eFactor = (($absDiff*$absDiff)/$bandWidth)/-2;
						
						$eFigure = $e**$eFactor;
						$gaussFigure = $oneOverRootTwoPi*$eFigure;
						
						$fracFactor = $absDiff/$bandWidth;
						$interTotal = $fracFactor*$gaussFigure;
											
						$probTotal = $probTotal + $interTotal;
					}
			}
		
		$actualTotal = $probTotal/($nP*$bandWidth);
		$actualTotal = sprintf("%.12f", $actualTotal);

		$actualTotal = "0.000000000001" if $actualTotal < 0.000000000001;

		return $actualTotal;

	}


# sub peakCaller takes two arguments: the name of the file to call peaks on, and the threshold density to call peaks at

sub peakCall
	{
		my $filename = $_[0];
		my $thresholdPeakCall = $_[1];
		my $centralBinSize = $_[2];
		my $binSize = $_[3];
		my $genomeSizeFile =  $_[4];
		my $stepSize =  $_[5];
		my %pValueTable = %{$_[6]};
		my $maxPeakPValue = $_[7];
	
		die "threshold calculation failed in peakCall sub.\n\n" if $thresholdPeakCall == 0;

		my %reads = ();
		my @line = ();
		my $token = 0;
		my $count = 0;
		my $binTotal = 0;
		my $backgroundTotal = 0;

		my @headerCheck = ();
		my $headerLineCheckerOff = 0;
		my @currentLineSplit = ();
		my $chromosome = 0;
		my $chromosomeSize = 0;
		my $currentBinStart = 0;
		my $currentBinEnd = 0;

		my $centralBinStart = 0;
		my $centralBinEnd = 0;
		my $centralBinScore = 0;
		my $binScore = 0;

		my $startString = 1;
		my @readStartPositionArray = ();
		my $arraySize = 0;

		my $peakOn = 0;
		my $backgroundScore = 0;
		my $segmentScore = 0;
		my $peakChr = 0;
		my $peakStart = 0;
		my $peakEnd = 0;
		my $peakLine = 0;
		my $peakMax = 0;
		my $peakNumber = 0;
		my @peakArray = ();
		my $peakChromSize = 0;

		my $i = 0;
#find the highest value of key in pValueTableArray:
		my $maxPValueTableValue = 0;
		if ($maxPeakPValue eq "ON")
			{
				foreach my $tekkers (keys %pValueTable)
					{
						$maxPValueTableValue = $tekkers if $tekkers > $maxPValueTableValue;
					}
			}

		# create hash with chromkey and chromsize value
				my %hashChromSizes = ();
				open (CHROMSIZE, "<", "$genomeSizeFile"); #$genomeSizeFile");
				while (<CHROMSIZE>)
					{
						chomp;
						my ($keyChromName, $valChromSize) = split /\t/;
						$hashChromSizes{$keyChromName} = $valChromSize; 
					}
				close CHROMSIZE;
				
				
				### print wig header in new wig file
				print STDERR "\n";
				###start of density track making
				open (ORIGINALSAM, "<$filename");
				while (<ORIGINALSAM>)
					{
						if ($headerLineCheckerOff == 0)
							{
								@headerCheck = split("", $_);
								$headerLineCheckerOff = 1 if $headerCheck[0] ne "@";
								next if $headerLineCheckerOff == 0;
							}
						# from HERE
						next unless ($_ =~ /\b(chr\w*)\b\t\b(\w*)\b/);
						my $thisChrom = $1;
						my $thisPos = $2;
						next unless exists $hashChromSizes{$thisChrom};


						#@currentLineSplit = split(" ", $_);
													
						if ("$thisChrom" ne "$chromosome")
							{			
								$chromosome = "$thisChrom";
								#skip read it chromosome is not included in chromosome size file
								next unless exists $hashChromSizes{$chromosome};
								$chromosomeSize = $hashChromSizes{$chromosome};
								print STDERR "\tfinding peaks in $chromosome...\n";

								$currentBinStart = $startString; #+ $offSetSubDR;
								$currentBinEnd = $currentBinStart + $binSize;
								$centralBinStart = ($currentBinEnd/2 - $centralBinSize/2) - 0.5; 
								$centralBinEnd = $centralBinStart + $centralBinSize;
								@readStartPositionArray = ();
								
							}
					
							#what the densities are calculated
							push(@readStartPositionArray, $thisPos);
							while (("$thisPos" > "$currentBinEnd") && ($currentBinEnd < $chromosomeSize ))
								{	
									while ($readStartPositionArray[0]<$currentBinStart)
										{
											shift(@readStartPositionArray);
										} 
							#to HERE
									$arraySize = @readStartPositionArray;
									for ($i = 0; $i < "$arraySize"; $i += 1)
										{
											if (($readStartPositionArray[$i] >= $currentBinStart) && ($readStartPositionArray[$i] < $currentBinEnd) && ($currentBinStart < $currentBinEnd))
												{
														$binScore += 1;
												}
										
											if (($readStartPositionArray[$i] >= $centralBinStart) && ($readStartPositionArray[$i] < $centralBinEnd) && ($centralBinStart < $centralBinEnd))
												{
														$centralBinScore += 1;
												}

										}
											
									$backgroundScore = ($binScore/$binSize)*$centralBinSize; #can int or round this<---
									$backgroundScore = sprintf("%.0f", $backgroundScore);
									$backgroundScore = int($backgroundScore);
									$segmentScore = $centralBinScore - $backgroundScore;

									if (($segmentScore >= $thresholdPeakCall) && ($peakOn ==0))
										{
											$peakOn = 1;
											$peakStart = $centralBinStart;
											$peakMax = $segmentScore;
											$peakEnd = $centralBinEnd;
											$peakChr = $chromosome;
											$peakChromSize = $chromosomeSize; # HERE
										}		

									if (($segmentScore >= $thresholdPeakCall) && ($peakOn ==1))
										{
											$peakMax = $segmentScore if $segmentScore > $peakMax;
											$peakEnd = $centralBinEnd;
										}

									if (($segmentScore < $thresholdPeakCall) && ($peakOn == 1))
										{
											$peakNumber += 1;
											if ($maxPeakPValue eq "ON")
												{
													$peakMax = "$pValueTable{$peakMax}" if exists $pValueTable{$peakMax};
													$peakMax = "0.000000000001" if $peakMax > $maxPValueTableValue;
												}

											$peakEnd = $peakChromSize if $peakChr ne $chromosome; # HERE

											$peakLine = "$peakChr\t$peakStart\t$peakEnd\tpk$peakNumber\t$peakMax\t+\n";
											print STDOUT "$peakLine";
											$peakOn = 0;
										}
									
									if (($segmentScore < $thresholdPeakCall) && ($peakOn == 0))
										{

										}

									$currentBinStart += "$stepSize";
									$currentBinEnd += "$stepSize";
									$centralBinStart += "$stepSize";
									$centralBinEnd += "$stepSize";
									$binScore = 0;
									$centralBinScore = 0;
								}
					}
				close ORIGINALSAM;
				print STDERR "\n\t$peakNumber peaks found in $filename.\n\n";
				return $peakNumber;
				
	}


##### readsPerBase calculates the number of reads in a sam file, and the number of base positions covered. it takes one parameter: a sam file name. it only counts mapped reads (where a base position is given.
# sub readsPerBase
# 	{
# 		my $rpbFile = $_[0];
# 		my $rpbreadCount = 0;
# 		my $rpbbaseCount = 0;
# 		my $rpbHeaderOff = 0;
# 		my $rpbcurrentBasePosition = 0;

# 		my @rpbsplitArray = ();
# 		my @rpbheaderCheckerArray = ();
# 		print STDERR "calculating rpb score for $rpbFile...\n";
# 		open (FILE, "<$rpbFile");
# 		while (<FILE>)
# 			{
# 				if ($rpbHeaderOff == 0)
# 					{
# 						@rpbheaderCheckerArray = split("", $_);
# 						$rpbHeaderOff = 1 if $rpbheaderCheckerArray[0] eq "@";
# 						next if $rpbHeaderOff == 0;
# 					}
			
# 				@rpbsplitArray = split(" ", $_);
# 				$rpbbaseCount += 1 if $rpbsplitArray[3] != $rpbcurrentBasePosition;
# 				$rpbreadCount += 1;
# 				$rpbcurrentBasePosition = $rpbsplitArray[3];
# 			}
# 		close FILE;
# 		my $readsPerBaseRpb = $rpbreadCount/$rpbbaseCount;
# 		my $roundedReadsPerBase = sprintf("%.3f", $readsPerBaseRpb);
# 		print STDERR "reads per base score: $roundedReadsPerBase\n\n";
# 		return 	$readsPerBaseRpb;	
# 	}


sub randomReadSelection
	{
		die "no parameters sent to Random Read Selection subroutine: requires filename and number of reads required\n\n" unless @_;
		my $file = $_[0];
		my $absoluteReadNumber = $_[1];
		
		my $counter = 0;
		my $numerator = 0;
		my $denominator = 0;
		my $percentage = 0;
		my $random = 0;
		my $check = 0;
		my $headerOff = 0;

		my @headerCheckerArray = ();

		print STDERR "\tcounting reads...\n";
		
		open (FILE, "<$file");
		while (<FILE>)
			{
				if ($headerOff == 0)
					{
						@headerCheckerArray = split("", $_);
						$headerOff = 1 if $headerCheckerArray[0] ne '@';
					}
				
				if ($headerOff ==1)
					{
						$counter += 1;
					}
			}
		close FILE;	
		
		$headerOff = 0;
		$percentage = ($absoluteReadNumber/$counter)*100;
		### more zeros gives more precision in estimation of numbers of reads to select.
		$numerator = int($percentage*10000);
		$denominator = 1000000;
		
		my $roundedPercent = sprintf("%.1f", $percentage);
		
		print STDERR "\tselecting ~$roundedPercent% of $counter reads...\n";
		open (FILE, "<$file");
		while (<FILE>)
			{
				if ($headerOff == 0)
					{
						@headerCheckerArray = split("", $_);
						print STDOUT "$_" if $headerCheckerArray[0] eq '@';
						$headerOff = 1 if $headerCheckerArray[0] ne '@';
					}
				
				if ($headerOff ==1)
					{
						$random = int(rand(1000000));
						if ($random <= $numerator)
							{
								print STDOUT "$_";
								$check += 1;
							}
					}
			}
			
		close FILE;
		print STDERR "\t$check reads randomly selected.\n\n";
	}
	
sub samFilter
	{
		### ****simplified from standalone
		#stamp package: for header detail
				
		# ** HCS denotes ordering as per $genomeSizeFile, and requires $genomeSizeFile to function
		###	
		# access $genomeSizeFile to get new order for sorted file
		# write comment line to start of new sorted sam file
		# option to INCLUDE header from alternate sam file
		# if sam header is present, write it to new sorted file
		# if sam header absent, or present, add new @SQ lines
				
		#	takes a mapped unsorted sam file that has been sorted by chromosome position with unix command:
		#
		#	sort -n --key=4,5 filename.sam > sortedFilename.sam
		#	
		#		and prints the header from the unsorted filename to the new filtered sam file, but orders the chromosomes as they appear in $genomeSizeFile
		
		#		new @SQ tags are added (in addition to original stampy @SQ tags) to denote the new order of chromosomes.		
		
		my $headerFile = $_[0];
		my $sortedFile = $_[1];
		my $mapq = $_[2];
		my $uq = $_[3];
		
		my $currentTempFile = 0;
		my $headerOff = 0;
		my @headerCheckerArray = ();	
		my @splitArray = ();
		my @currentUQArray = ();		
		
		my $passedFilter = 0;
		my $totalReads = 0;
		my $i = 0;
		my $percentagePassedFilter = 0;
		my $numberOfChromosomes = 0;
		my $tempSwitch = 0;
		my $unmappedCount = 0;
		my $percentageUnmapped = 0;
		
		my $chromCount = 0;
		my @chromNames = ();
		my @chromSizes = ();
		my @line = ();
		my %chromReadCount = ();
		my $headOne = 0;
		my $samHeader = "FALSE";
		my $percentageDuplicates = 0;
		my $testString = 0;
		my $oldTestString = 0;
		my $duplicateCount = 0;
		
				
		print STDOUT "\@CO\ttime:$timeString\tprog:$0\tuser:$user\tcwd:$dir\tmapq:$mapq\tuq:$uq\n"; #inserts new header comment line with current additional details
		open(HEADERFILE, "<$headerFile");
		while (<HEADERFILE>) # then prints remainder of original header and stops at with first non-header line
			{
				@splitArray = split("", $_);
				if ($splitArray[0] eq '@')
					{
						$samHeader = "TRUE";
						print STDOUT "$_";
					}
				else {last;}
					
			}
		close HEADERFILE;
		
		print STDERR "\tsam header detected, original header transferred to new file. \n" if $samHeader eq "TRUE"; #confirms status of header
		print STDERR "\twarning: sam header not detected.\n" if $samHeader eq "FALSE";
		
		# create new @SQ tags in the order of $genomeSizeFile. These tags are indicated by @SQ\tSNN:<chrom.name>:hg19ChromSize\tLN:<chrom.size>. The ultrafiltered sam file will be made in this order
		print STDERR "\tcreating chromosome order list from $genomeSizeFile...\n";
		open (CHROMSIZES, "<$genomeSizeFile");
		while (<CHROMSIZES>)
			{
				@line = split(" ", "$_");
				push (@chromNames, "$line[0]");
				push (@chromSizes,  "$line[1]");
			}
			
		# initialise chromosome passed filter counts hash with zero value for each:
		for ($i = 0; $i < @chromNames; $i += 1)
			{
				$chromReadCount{$chromNames["$i"]} = 0;
			}
		
		$currentTempFile = "$headerFile" . ".temp";
		
		# makes a new temporary file containing only reads that pass the mapq and uq filters
		print STDERR "\tundertaking initial mapq and uq filtering...\n";
		open (SORTEDINPUT, "<$sortedFile");
		open (TEMP, ">$currentTempFile");
		
		my %presentChroms = ();
		my @chromsToPrint = (); 
		my @testChroms = ();
		
		while (<SORTEDINPUT>)
			{
				if ($headerOff == 0)
					{
						@headerCheckerArray = split("", $_);
						$headerOff = 1 if $headerCheckerArray[0] ne '@';
					}
				
				if ($headerOff ==1)
					{
						@splitArray = split(" ", $_);
						if ($pcrDuplicates eq "ON") # next seven lines added 20130507 to put in PCR duplicate checker
							{
								$oldTestString = $testString;
								$testString = "$splitArray[2]"."$splitArray[3]"."$splitArray[4]"."$splitArray[8]";
								$duplicateCount += 1 if $testString eq $oldTestString;
								next if $testString eq $oldTestString;
							}
						
						$unmappedCount += 1 if $splitArray[3] == 0;
						if ($splitArray[3] > 0)
							{
								@currentUQArray = split(":", $splitArray[11]);
								if (($splitArray[4] >= $mapq) && ($currentUQArray[2] <= $uq))  
									{
										$presentChroms{"$splitArray[2]"} = 0 unless exists $presentChroms{"$splitArray[2]"}; #makes a hash containing as keys all chromosome names presnt
										print TEMP "$_";
										$passedFilter += 1;
										$chromReadCount{$splitArray[2]} += 1;
									}
							}
						$totalReads += 1;
					}		
			}
		close SORTEDINPUT;
		close TEMP;
		
		@testChroms = %presentChroms; # makes an array of chromosomes that have been found in the file
		for ($i = 0; $i < @chromNames; $i += 1) # adds to a new array chromosomes found in the file, in the order they appear in the chrom.sizes file
			{
				push (@chromsToPrint, "$chromNames[$i]") if exists $presentChroms{"$chromNames[$i]"};
			}
				
		# adding @SQ lines to sorted sam header that contain number of reads that have passed the filters	
		my $accuratePercent = $passedFilter/$totalReads*100;
		$accuratePercent = sprintf("%.2f", $accuratePercent);
		
		for ($chromCount = 0; $chromCount < @chromNames; $chromCount += 1)
			{
				print STDOUT "\@SQ\tSN:$chromNames[$chromCount]:readsPassedFilter:$chromReadCount{$chromNames[$chromCount]}\tLN:$chromSizes[$chromCount]\n";
			}
		
		print STDOUT "\@CO\ttotalReads:$totalReads\treadsPassedFilter:$passedFilter\t percentPassed:$accuratePercent%\tmapQpoint:$mapq\tUQpoint:$uq\n";
		print STDERR "\tnew header completed.\n\tbeginning ordering of chromosomes...\n\n";
		
		## this is the ordering section		
		$numberOfChromosomes = @chromsToPrint;
		for ($i = 0; $i < $numberOfChromosomes; $i += 1)
			{
					print STDERR "\tsorting $chromsToPrint[$i]...\n";
					open (TEMPA, "<$currentTempFile");
					while (<TEMPA>)
							{
								print STDOUT "$_" if ($_ =~ m#\b$chromsToPrint[$i]\b#);
							}
	
					close TEMPA;
			}
		
		unlink $currentTempFile;
				
		$percentagePassedFilter = ($passedFilter/$totalReads)*100;
		$percentageUnmapped = ($unmappedCount/$totalReads)*100;
		$percentagePassedFilter = sprintf("%.2f", $percentagePassedFilter);
		$percentageUnmapped = sprintf("%.2f", $percentageUnmapped);
		$percentageDuplicates = ($duplicateCount/$totalReads)*100;
		$percentageDuplicates = sprintf("%.2f", $percentageDuplicates);
		print STDERR "\n\t$percentageUnmapped% of total reads were unmapped.\n";
		print STDERR "\t$duplicateCount PCR duplicates ($percentageDuplicates% of total) were removed.\n" if $pcrDuplicates eq "ON";
		print STDERR "\t$passedFilter ($percentagePassedFilter%) reads of total $totalReads successfully passed filter.\n\n";
	}	


sub densityAnalyzer
	{
		
		#removed requirement for HCS passed filter
		my $filename = $_[0];
		my $stepSizeSubDR = $_[1];
		my $tailSubDR = $_[2];
		our $sigmaSubDR = $_[3];
		my $threshSubDR = $_[4];
		my $offSetSubDR = $_[5];
		my $maxThreshSubDR = $_[6];
		my $genomeSizeFile = $_[7];

		my $bandWidth = $tailSubDR*2 +1;
		my $binSize = $tailSubDR*2 +1;
		my $span = $stepSizeSubDR;
		my $gaussCount = 0;
		my $headerLineCheckerOff = 0;
		my $currentBinStart = 0;
		my $currentBinEnd = 0;
		my $positionCount = 0;
		my $i = 0;
		my $offStart;
		my $offEnd;
		my $chromosome = 0;
		my $startString = 1;
		my $binScore = 0;
		my $arraySize = 0;
		my $chromosomeSize = 0;
		my $headerSwitch = 0;
		my $gaussOutput = 0;
		my $memSaverGauss = 0;
		my $sum = 0;
		my $roundedFigure = 0;
		my $gapCount = 0;
		my $stub = $filename;chomp $stub;for ($i = 0; $i < 4; $i += 1){chop $stub;}
		my $trackName = "$stub.step$stepSizeSubDR.tail$tailSubDR.sig$sigmaSubDR.thresh$threshSubDR.off$offSetSubDR";
		my $trackDescription = "$stub" . "-$stepSizeSubDR.$tailSubDR.$sigmaSubDR.$threshSubDR.$offSetSubDR" . "-";
		
		my @currentLineSplit = ();
		my @readStartPositionArray = ();
		my @headerCheck = ();
		my @splitArray = ();
		my @gaussArray = ();
		our @fixedGaussArray = ();
		my @gaussSubArray = ();
		
		# generates the gaussian values based on sigma and tailsize; eg for tail of 150, 301 values are generated and placed in the fixedGauss reference array
		for ($i = -$tailSubDR; $i <= $tailSubDR; $i += 1)
				{
					$gaussOutput = $sigmaSubDR*gaussEngine($i);
					push(@fixedGaussArray, $gaussOutput);
				}
		
		# if thresh is not set by user, calculates default bin thresh # HERE (moved lines 1750-1758 up)
		my %hashChromSizes = ();
		open (CHROMSIZE, "<$genomeSizeFile");
		while (<CHROMSIZE>)
			{
				chomp;
				my ($keyChromName, $valChromSize) = split /\t/;
				$hashChromSizes{$keyChromName} = $valChromSize; 
			}
		close CHROMSIZE;

		# if ($threshSubDR eq "DEFAULT") # HERE
		# 	{	# HERE
		my $threshCalcHeaderOff = 0; #<-removed == from here 20130703
		my $threshCalcReadCount = 0;
		my @threshCalcHeader = ();
		my $availableBaseCount = 0;
		print STDERR "\tcalculating default minimum bin threshold...\n"; # HERE
		$genomeBaseCount = baseCountGenomeSizeFile($genomeSizeFile);
		my %chromsPresent;
		open (THRESHCALC, "<$filename");
		while (<THRESHCALC>) 
			{
				if ($threshCalcHeaderOff == 0)
					{
						@threshCalcHeader = split("", $_);
						$threshCalcHeaderOff = 1 if $threshCalcHeader[0] ne "@";
						next if $threshCalcHeaderOff == 0;
					}
				my $lineChrom = $1 if $_ =~ /(\bchr\w*\b)/; # HERE
				$availableBaseCount = $availableBaseCount + $hashChromSizes{$lineChrom} unless defined $chromsPresent{$lineChrom};
				$chromsPresent{$lineChrom} = 1 if defined $hashChromSizes{$lineChrom}; # HERE
				#print STDERR "$lineChrom\n" if defined $hashChromSizes{$lineChrom};  # HERE
				$threshCalcReadCount += 1 if defined $hashChromSizes{$lineChrom}; # HERE
			}
		die "\tnone of the chromosomes listed in $genomeSizeFile have been detected in $filename.\n\n" unless scalar keys %chromsPresent >= 1;
		my $threshSetting = $threshSubDR; # HERE
		my $threshCalced = 2*(($threshCalcReadCount/$availableBaseCount)*$binSize);# HERE
		$threshCalced = sprintf("%.0f", $threshCalced); # HERE
		$threshSubDR = $threshCalced if $threshCalced eq "DEFAULT"; # HERE
		$threshSubDR = sprintf("%.0f", $threshSubDR);
		print STDERR "\tgenome size file: $genomeSizeFile\n\tgenome total size: $genomeBaseCount\n\tavailable bases: $availableBaseCount\n\ttotal number of reads: $threshCalcReadCount\n"; # HERE
		print STDERR "\tdefault calculated minimum bin threshold: $threshCalced\n\n"; # HERE
		print STDERR "\tusing assigned threshold: $threshSubDR\n\n" if $threshSetting ne "DEFAULT"; # HERE
			#} # HERE
			
		# create hash with chromkey and chromsize value
		
		
		### print wig header in new wig file
		print STDOUT "track type=wiggle_0 name=$trackName description=$trackDescription viewLimits=0:500 autoScale=off gridDefault=on color=0,160,255 maxHeightPixels=60:60:11 visibility=full windowingFunction=mean+whiskers\n";
		
		###start of density track making
		open (ORIGINALSAM, "<$filename");
		while (<ORIGINALSAM>)
			{
				if ($headerLineCheckerOff == 0)
					{
						@headerCheck = split("", $_);
						$headerLineCheckerOff = 1 if $headerCheck[0] ne "@";
						next if $headerLineCheckerOff == 0;
					}

				#from HERE
				next unless ($_ =~ /\b(chr\w*)\b\t\b(\w*)\b/);
						my $thisChrom = $1;
						my $thisPos = $2;
						next unless exists $hashChromSizes{$thisChrom};
				#@currentLineSplit = split(" ", $_);
											
				if ("$thisChrom" ne "$chromosome")
					{			
						$chromosome = "$thisChrom";
						#next unless exists $hashChromSizes{$chromosome};
						$chromosomeSize = $hashChromSizes{$chromosome};
											
						print STDERR "\tmaking density track for $chromosome...\n"; # HERE
						
						$currentBinStart = $startString; #+ $offSetSubDR;
						$currentBinEnd = $currentBinStart + $binSize;
						@readStartPositionArray = ();
						@gaussArray = ();
						$positionCount = 0;
						$gapCount = 1;
					}
			
					#where the peak calling is decided:
					push(@readStartPositionArray, $thisPos);
					$offEnd = $currentBinEnd + $offSetSubDR;
					while (("$thisPos" > "$currentBinEnd") && ($offEnd < $chromosomeSize ))
						{	
							while ($readStartPositionArray[0]<$currentBinStart)
								{
									shift(@readStartPositionArray);
								} 
			#to HERE
							$arraySize = @readStartPositionArray;
							for ($i = 0; $i < "$arraySize"; $i += 1)
								{
									if (($readStartPositionArray[$i] >= $currentBinStart) && ($readStartPositionArray[$i] < $currentBinEnd) && ($currentBinStart < $currentBinEnd))
										{
												$binScore += 1;
										}
								}
											
							if (($binScore > $threshSubDR) && ($binScore <= $maxThreshSubDR))
								{
									
									for ($gaussCount = 0; $gaussCount < $bandWidth; $gaussCount += 1)
										{
											$gaussArray[$gaussCount] = 0 if !defined $gaussArray[$gaussCount];
											$gaussArray[$gaussCount] = $gaussArray[$gaussCount] + $binScore*$fixedGaussArray[$gaussCount] if ((defined $gaussArray[$gaussCount]) && (defined $fixedGaussArray[$gaussCount]));
											#print STDOUT "$gaussArray[$gaussCount]:$binScore:$fixedGaussArray[$gaussCount]\n";
										}
								}
												
							@gaussSubArray = splice(@gaussArray, 0, "$stepSizeSubDR");
							
							$sum += $_ for (@gaussSubArray);
							$roundedFigure = int(10*$sum);
							
							if ($roundedFigure == 0)
								{
									$gapCount = 1;
									@gaussArray = ();
								}
							
							if (($roundedFigure >0) && ($gapCount == 0))
								{
									for ($gaussCount = 0; $gaussCount < $stepSizeSubDR; $gaussCount += 1)
										{
											$memSaverGauss = int(10*$gaussSubArray[$gaussCount]);
											print STDOUT "$memSaverGauss\n" if $offStart > 0;
										}
								}
							
							if (($roundedFigure > 0) && ($gapCount > 0))
								{
									$offStart = $currentBinStart + $offSetSubDR;
									if ($offStart >0)
										{
											print STDOUT "fixedStep chrom=$chromosome start=$offStart step=1\n";
									
											for ($gaussCount = 0; $gaussCount < $stepSizeSubDR; $gaussCount += 1)
												{
													$memSaverGauss = int(10*$gaussSubArray[$gaussCount]);
													print STDOUT "$memSaverGauss\n";
												}
										}
										
									$gapCount = 0;
								}
							
							$sum = 0;			
							$positionCount += 1;
							$currentBinStart += "$stepSizeSubDR";
							$currentBinEnd += "$stepSizeSubDR";
							$binScore = 0;                 
						}
			}
		close ORIGINALSAM;
		print STDERR "\n\tDensity track created.\n\n";
	}


sub peakSignalIntensity
	{
		my $file = $_[0];
		my $distanceFromPeakCentre = $_[1];
		my $numberOfPeaks = $_[2];
		my $targetSamFile = $_[3];
		
		my $counter = 0;
		my $numerator = 0;
		my $denominator = 0;
		my $percentage = 0;
		my $random = 0;
		my $check = 0;
		my $actualPeakCount = 0;
		my $reCalc = 1;
		my $currentTargetPeak = 0;		
		my $lowPeak = 0;	
		my $highPeak = 0;		
		my $peakCentre = 0;
		my $peakChrom = 0;
		my $peakName = 0;
		my $lowLimit = 0; 
		my $highLimit = 0;
		my $readValue = 0;
		my $peaksCompleted = 0;
		my $currentValue = 0;
		my $currentXCoord = 0;
		my $totalPositions = 0;
		my $headerLineCheckerOff = 0;
		my $s = 0;

		my @meanGlobalBin = ();
		my @globalBin = ();
		my @splitArray = ();
		my @headerCheck = ();
		my @bedPeaksArray = ();
		my @splitSamArray = ();

		#open target bed file, and randomly select $absoluteReadNumber reads, and put these into an array new file

		open (FILE, "<$file");
		while (<FILE>)
			{
					$counter += 1;
			}
		close FILE;	
		
		$percentage = ($numberOfPeaks/$counter)*100;
		$numerator = int($percentage*100000);
		$denominator = 10000000;
		my $roundedPercent = sprintf("%.2f", $percentage);
		
		open (FILE, "<$file");
		while (<FILE>)
			{
				@splitArray = split(" ", $_);
					$random = int(rand(10000000));
					if ($random < $numerator)
						{
							chomp $_;
							$bedPeaksArray[$actualPeakCount] = $_;
							$actualPeakCount += 1;
						}
			}
		close FILE;
		
		print STDERR "\ttotal peaks: $counter; percentage for $numberOfPeaks peaks: $roundedPercent%\n\t$actualPeakCount peaks were randomly selected.\n";
		print STDERR "\tcompleted target peak survey.\n\tcalculating average signal intensity for $actualPeakCount peaks...\n\n";
			
		open (SAMFILE, "<$targetSamFile");
		while (<SAMFILE>)
			{
				if ($headerLineCheckerOff == 0)
					{
						@headerCheck = split("", $_);
						next if $headerCheck[0] eq "@";
						$headerLineCheckerOff = 1 if $headerCheck[0] ne "@";
					}
				
				if ($headerLineCheckerOff == 1)
					{
						@splitSamArray = split("\t", $_);
												
						if ($reCalc == 1)
							{
								$currentTargetPeak = $bedPeaksArray[$currentValue];
								@splitArray = split("\t", $currentTargetPeak);
								$lowPeak = $splitArray[6];
								$highPeak = $splitArray[7];
								$peakCentre = int(($lowPeak + $highPeak)/2);
								$peakChrom =  $splitArray[0];
								$peakName =  $splitArray[3];
								
								$lowLimit = $peakCentre - $distanceFromPeakCentre;
								$highLimit = $peakCentre + $distanceFromPeakCentre;
								
								last if $lowPeak == 0;
								$reCalc = 0;
							}
										
						if (($splitSamArray[2] eq "$peakChrom") && ($splitSamArray[3] >= "$lowLimit") && ($splitSamArray[3] <= "$highLimit"))
							{
									$readValue = $splitSamArray[3] - $lowLimit;
									$globalBin[$readValue] += 1;
							}

						if ((($splitSamArray[2] eq "$peakChrom") && ($splitSamArray[3] > "$highLimit")) || (($splitSamArray[2] ne "$peakChrom") && ($splitSamArray[3] < "$lowLimit")))
							{
								$reCalc = 1;
								$currentValue += 1;
								$peaksCompleted += 1;
							}
					}
			}
						
		close SAMFILE;

		$currentXCoord = -1*$distanceFromPeakCentre;
		$totalPositions = ($distanceFromPeakCentre*2) + 1;
		for ($s = 0; $s < $totalPositions; $s += 1)
				{
					$meanGlobalBin[$s] = $globalBin[$s]/$peaksCompleted;
					print STDOUT "$currentXCoord\t$meanGlobalBin[$s]\n";
					$currentXCoord += 1;
				}

		print STDERR "\t$peaksCompleted peaks completed.\n\n";
	}

sub topPeaks
	{
		
		my $file = $_[0];
		my $numberToPrint = $_[1];
		my @line = ();
		my %largeSORTHash = ();
		my $tempCounter = 0;
		my @orderedArray;
		my $counter = 0;
		my @individualArray;
		my @totalArray;
		my $x;
		my $limit;
		
		print STDERR "\tsorting $numberToPrint peaks in $file...\n\n";

		
		open (HASHSORT, "<", "$file") || die  "Unable to open $file\n\n";;
		while (<HASHSORT>)
				{
					@line = split(" ", $_);
					push (@{$largeSORTHash{$line[4]}}, $_);
				}
		close HASHSORT;
		
		foreach my $key (reverse sort {$a <=> $b} keys %largeSORTHash) 
			{
				push (@orderedArray,  "@{$largeSORTHash{$key}}");
				$counter += 1;
			}
		
		for (@orderedArray)
			{
				@individualArray = split("\n", $_);
				for (@individualArray)
					{
						$_ =~ s/^\s//; 
						push (@totalArray, $_);
					}
			}
		my $count = @totalArray;
		$limit = $numberToPrint;
		$limit = $count if $numberToPrint eq 'ALL';
		for ($x = 0; $x < $limit; $x += 1)
			{
				print STDOUT "$totalArray[$x]\n";
			}

		
		print STDERR "\tsorted.\n\n";
		return;

# old method involves temp files. new is better.
		# my $file = $_[0];
		# my $numberToPrint = $_[1];
		

		# my @line = ();
		# my %largeSORTHash = ();
		# my $tempList = "$file" . ".temp";
		# my $tempCounter = 0;
		
		# print STDERR "sorting $numberToPrint peaks in $file...\n\n";
			
		# open (HASHSORT, "<$file");
		# while (<HASHSORT>)
		# 		{
		# 			@line = split(" ", $_);
		# 			push @{$largeSORTHash{$line[4]}}, $_;
		# 			#$largeSORTHash{$line[4]} = "$_";
		# 		}

		# open (TEMPLIST, ">$tempList");
		# foreach my $key (reverse sort {$a <=> $b} keys %largeSORTHash) 
		# 	{
		# 		print TEMPLIST @{$largeSORTHash{"$key"}};
		# 	}
		# close TEMPLIST;
		
		# open (TEMPLIST, "<$tempList");
		# while (<TEMPLIST>)
		# 	{
		# 		print STDOUT "$_";
		# 		$tempCounter += 1;
		# 		unless ($numberToPrint eq "ALL") 
		# 			{
		# 				last if $tempCounter == $numberToPrint;
		# 			}
		# 	}
		# close TEMPLIST;
		# unlink $tempList;
	}

sub samNSort
	{
		my $file = $_[0];
			
		my @line = ();
		my %largeSORTHash = ();
		my $currentKey = 0;
		my $maxValue = 0;
		my $headerOff = 0;
		my @headerCheckerArray;
		my $count;		
		print STDERR "\tsorting $file...\n\n";
		 
		open (HASHSORT, "<$file");
		while (<HASHSORT>)
				{
					if ($headerOff == 0)
						{
							@headerCheckerArray = split("", $_);
							print STDOUT "$_" if $headerCheckerArray[0] eq '@';
							$headerOff = 1 if $headerCheckerArray[0] ne '@';
						}
				
					if ($headerOff ==1)
						{
							@line = split(" ", $_);
							push @{$largeSORTHash{$line[3]}}, $_;
							$maxValue = $line[3] if $line[3] > $maxValue;
							$count += 1;
							if ($count%1000000 == 0)

								{
									print STDERR "\tcollected $count reads...\n";
									
								}
						}
				}
		my $number = keys %largeSORTHash;
		
		for ($currentKey = 0; $currentKey <= $maxValue; $currentKey += 1)
			{
				print STDOUT @{$largeSORTHash{"$currentKey"}} if defined $largeSORTHash{"$currentKey"};
			}

		print STDERR "\tsorted.\n\n";
	}
	
sub gaussEngine
	{		
		my $in = "@_";
		my $denominator = 1/sqrt(2*3.14159265*$sigma**2);
		my $sigmaSubGEsquaredTwice = ($sigma*$sigma)*2;
		my $power = -($in**2)/$sigmaSubGEsquaredTwice;
		my $eToPower = 2.7182818284590451**$power;
		my $gaussed= $eToPower*$denominator;
		my $roundGaussed = sprintf("%.5f", $gaussed);
		return $roundGaussed;
	}

sub baseCountGenomeSizeFile
	{	
		my $gSizeFile = $_[0];
		my $gTotalSize = 0;
		my @gline = ();
		open (GSIZE, "<$gSizeFile");
		
		while (<GSIZE>)
			{
				@gline = split(" ", $_);
				$gTotalSize += $gline[1];
			}
		return $gTotalSize;
	}

sub fileTest
	{
		my $format = $_[0];
		my $file = $_[1];

		
		my $headerOff = 0;
		my @line;
		my @sameLine;
		my $samFORMAT;
		my $samDisorder;
		my $samOrder = 0;
		my $header;
		my $samCount = 0;
		my $chrOrder;
		my $chrDisorder;
		my %samTags = ('@CO' => 1, '@RG' => 1, '@SQ' => 1, '@HD' => 1, '@PG' => 1);

		### defining SAM format
		if ($format eq 'SAM')
			{
				$header = 'NONE';
				$samFORMAT = 'OK';
				$samDisorder = 'OK';
				$chrDisorder = 'OK';
				
				if (-z $file)
					{
						$samFORMAT = 'BAD';
						$samDisorder = 'BAD';
						$header = 'NONE';
						$chrDisorder = 'BAD';
					}
				
				open (my $samFile, "<", $file) or die "Unable to open $file.\n\n";
				while (<$samFile>)

					{
						next unless $_ =~ /\w/;
						if (eof)
							{
								$samFORMAT = 'BAD';
								$samDisorder = 'BAD';
								$header = 'NONE';
								$chrDisorder = 'BAD'
							}

						
						if ($headerOff == 0)
							{
								@line = split(" ", $_);
								$header = 'OK' if defined $samTags{$line[0]};
								$headerOff = 1 unless exists $samTags{$line[0]};
								#print STDERR "Header off = $headerOff\n";
								#$headerOff = 1 if !defined $line[0];
							}
						if ($headerOff == 1)
							{
								
								@sameLine = split(" ", $_);

								$samFORMAT = 'BAD' unless defined $sameLine[7];
								if (defined $sameLine[7])
									{
										if (
											($sameLine[0] =~ /[^!-?A-~]/) ||
											($sameLine[1] =~ /\D/) ||
											(($sameLine[2] !~ /chr/) && ($sameLine[2] ne '*')) ||
											($sameLine[3] =~ /\D/) ||
											($sameLine[4] =~ /\W/) ||
											(($sameLine[5] =~ /\W/)  && ($sameLine[2] ne '*')) ||
											($sameLine[6] =~ /[^!-?A-~]/) ||
											($sameLine[7] =~ /[^!-?A-~]/) ||
											($#sameLine < 8)
										)
											{
												$samFORMAT = 'BAD';
											}
									} 
									
								$samDisorder = 'BAD' unless defined $sameLine[2];
								$chrDisorder = 'BAD' unless defined $sameLine[3];

								if ((defined $sameLine[3]) && ($samCount >= 1))
									{
										if ($sameLine[3] =~ /\d/)
											{
												unless ($sameLine[3] >= $samOrder)
													{
														$samDisorder = 'BAD' unless $sameLine[2] ne $chrOrder;
													}
											}
										else
											{
												$samDisorder = 'BAD'
											}
									}
								
								if ((defined $sameLine[2]) && ($samCount >= 1))
									{
										unless ($sameLine[2] eq $chrOrder)
											{
												$chrDisorder = 'BAD';
											}
									}
								

								$samOrder = $sameLine[3] if defined $sameLine[3];
								$chrOrder = $sameLine[2] if defined $sameLine[2];
								$samCount += 1;
								last if $samCount == 20; 
							}
					}

				return ($samFORMAT, $samDisorder, $header, $chrDisorder);

			}

		if ($format eq 'BED')
			{
				
				my $bedFORMAT = 'OK';
				my $bedCount = 0;
				
				if (-z $file)
					{
						$bedFORMAT = 'BAD';
					}
				open (my $bedFile, "<", $file) or die "Unable to open $file\n\n";

				while (<$bedFile>)
					{
						@line = split (" ", $_);
						
						if (eof)
							{
								$bedFORMAT = 'BAD';
							}


						if (
								($line[0] !~ /chr/) ||
								($line[1] =~ /\D/) ||
								($line[2] =~ /\D/) 
							)
								{
									$bedFORMAT = 'BAD';
								}

						$bedCount += 1;
						last if $bedCount == 2;
					}

				return ($bedFORMAT);

			}

		if ($format eq 'CHROM')
			{
				
				my $chromFORMAT = 'OK';
				my $chromCount = 0;
				
				if (-z $file)
					{
						$chromFORMAT = 'BAD';
					}
				open (my $chromFile, "<", $file) or die "Unable to open $file\n\n";

				while (<$chromFile>)
					{
						@line = split (" ", $_);
						
						# HERE
						# if (eof)
						# 	{
						# 		$chromFORMAT = 'BAD';
						# 	}


						if (
								($line[0] =~ /[^!-?A-~]/) ||
								($line[1] =~ /\D/) ||
								($#line != 1)

							)
								{
									$chromFORMAT = 'BAD';
								}

						$chromCount += 1;
						last if $chromCount == 1;
					}

				return ($chromFORMAT);

			}
	}

sub usage
	{
		my $indentLevel;
		my $count;
		while (<DATA>)
			{
				if ($_ =~ m/^=head1/)
					{
						$_ =~ s/=head1 //;
						$indentLevel = '';
					}

				if ($_ =~ m/^=head2/)
					{
						$_ =~ s/=head2 //;
						$indentLevel = '   ';
					}

				if ($_ =~ m/^=head3/)
					{
						$_ =~ s/=head3 //;
						$indentLevel = '      ';
					}

				last if $_ =~ m/^=cut/;

				if ($_ =~ m/I</ .. />/)
					{
						
						$_ =~ s/I<//;
						$_ =~ s/>//;
					}

				print STDOUT "$indentLevel" . "$_";

			}

	}

__DATA__

=head1 NAME

	PeaKDEck.pl v1.1 - a kernel density estimator based peak caller for I<DNaseI-seq> data. 

=head1 SYNOPSIS

	Numerical sorting:
	PeaKDEck.pl –NS mappedData.sam > sortedData.sam

	Filtering: 
	PeaKDEck.pl –F sortedData.sam –g hg19.chrom.sizes –q 15 –PCR ON > filteredData.sam

	Random read selection: 
	PeaKDEck.pl –R filteredData.sam –nr 20000000 > filteredData.20m.sam

	Density analysis:
	PeaKDEck.pl –D filteredData.20m.sam –g hg19.chrom.sizes > densityTrack.wig

	Peak calling:
	PeaKDEck.pl –P filteredData.20m.sam –g hg19.chrom.sizes > peakList.bed

	Size ordering of peaks: 
	PeaKDEck.pl –T peakList.bed –n 50000 > orderedPeakList.bed

=head1 ARGUMENTS

=head2 Numerical sorting (-NS I<mappedReads.sam>)
	
	Sorts sam format reads by base start position, irrespective of chromosome. Equivalent
	to the linux/unix/osx command: [sort -n --key=4,5 filename.sam > sortedFilename.sam].
	For fast results, memory corresponding to ~2.5 times file size should be available.

=head2 Sam filtering (-F I<mappedReads.sam>)
	
	Sorts sam files by chromosome, in the order that chromosomes appear in the chromosome
	size file. The chromosome size file is mandatory. The chromsome size file is a plain 
	text, tab separated file in the format: 

		chr1	249250621
		chr2	243199373
		chr3	198022430
		.................
		chrN 	size(bp)

=head3 Mandatory settings

	-g /path/to/chromosomeSizeFile.txt
		Specifies the path to the text file containing tab separated list of chromosome names
		and sizes.

=head3 Optional settings

	-q integerValue
		Specifies a mapq cutoff score for filtering. Reads with a mapq score less than the 
		supplied value will be removed from the resulting filtered file. By default, -q is 
		set to zero, so no filtering for mapq scores will occur.

	-u integerValue
		Specifies a UQ base mismatch score for filterering. Reads with mismatch scores greater
		than this value will be removed from the filtered dataset. By default, -u is set to 
		10000, so that no filtering by uq score will occur.

	-i samHeaderFile.sam
		Specifies a file containing a sam header, which if set, will be included at the beginning
		of the newly filtered file.

	-PCR ON|OFF
		Allows PCR duplicate reads to be removed from sam file. Reads are considered PCR duplicates
		if adjacent reads have identical chromosome, start position, mapq score, and sequence. By 
		default -PCR is set to OFF, so no filtering of PCR duplicates will occur. To detect PCR
		duplicates, chromosomes must be in numerical order (see Numerical sorting above).

=head2 Random read selection (-R I<mappedReads.sam>)
	
	Randomly selects a target number of reads from a specified sam file. Selected reads are
	printed to STDOUT by default.

=head3 Mandatory settings

	-nr integer
		Specifies the target number of reads to be randomly selected from the given sam file.
		The number of reads must by a positive whole number.

=head2 Density analyzer (-D I<mappedOrderedReads.sam>)

	Creates a smoothed, unitless read density track in wig format, representing the distribution 
	of reads in the given sam file. Sam files must be grouped by chromosome, and ordered by read 
	start position (see Numerical sorting and Sam filtering above). The order of chromosomes in 
	the density track is determined by the order in which they appear in the mandatory 
	chromosome size file (see Sam filtering for chromosome size file format). By default, the 
	results are printed to STDOUT.


=head3 Mandatory settings

	-g /path/to/chromosomeSizeFile.txt
		Specifies the path to the text file containing tab separated list of chromosome names
		and sizes.

=head3 Optional settings

	-n positiveInteger
		Specifies the one-tailed size of the smoothing bin. By default, -t is set to 150, giving a bin
		size of 300 bp. This value determines both the size of sampling bin, and the width of the 
		Gaussian probability density function used to calculate read densities, and must by a 
		postitive whole number.

	-STEP positiveInteger
		Specifies the size of steps by which the probability density function and sampling bin move 
		along the genome. By default, -STEP is set to 100. Smaller step sizes proportionately
		increase the number of calculations carried out, and therefore the time taken for the 
		analysis. -STEP must be a positive whole number.

	-d positiveInteger
		Specifies the standard deviation of the probability density function. This value determines how 
		broadly the read density scores are spread over each sampling bin, and therefore determines
		the degree of smoothing that occurs. By default -d is set to 50, and must be a positive whole
		number.

	-t positiveInteger
		Specifies a low threshold, below which read density scores won't be included in probability
		density function calculations. By default, -t is set to the number of reads expected to occur
		in the set bin size if the number of reads in the dataset were randomly distributed. All reads
		present in the data set will be included in the analysis if -t is set to 0. -t must be a non-
		negative whole number.

	-m positiveInteger
		Specifies a high threshold, above which read density scores won't be included in probability
		density function calculations. By default, -m is set to 100000000, ensuring that no reads
		will be excluded from analysis in default settings.

	-o integer
		Specifies a track offset. All positions in the resulting wig file will be offset by this value.
		For DNaseI-seq data, the read start sites are considered DNaseI cutting sites, and so by 
		default, -o is set to 0. If the centre of the DNA fragment is considered the point of interest
		(for example, in ChIP-seq), setting -o to half the average fragment size may give a more
		precise depiction of signal localisation.

=head2 Peak calling (-P I<mappedOrderedReads.sam>)

	Identifies peaks in the provided sam file, and provides output in bed format to STDOUT. Sam 
	files must be grouped by chromosomes, and ordered by read position (see Numerical sorting 
	and Sam filtering above). The order of chromosomes in the peak file is determined by 
	the order in which they appear in the mandatory chromosome size file (see Sam filtering for 
	chromosome size file format).

=head3 Mandatory settings

	-g /path/to/chromosomeSizeFile.txt
		Specifies the path to the text file containing tab separated list of chromosome names
		and sizes.

=head3 Optional settings

	-bin positiveInteger
		Specifies the size of the central sampling bin. By default, -bin is set to 300, which 
		represents the expected average feature size. -bin must by set to a positive whole number

	-back positiveInteger
		Specifies the size of the background sampling bin. By default, -back is set to 3000, ten
		times the size of the central sampling bin. -back must be set to a positive whole number
		and must be larger than the size of the central bin.

	-STEP positiveInteger
		Specifies the size of steps by which the sampling bin moves along the genome. By default, 
		-STEP is set to 100. Smaller step sizes proportionately increase the number of calculations 
		carried out, and therefore the time taken for the analysis. -STEP must be a positive whole 
		number.

	-FLAT positiveInteger
		Specifies a flat threshold for peak calling in reads per bin. When -FLAT is set, the
		threshold calculated by PeaKDEck for peak calling is overridden, and the value given 
		by -FLAT is used in its place. FLAT must by a positive number.

	-b /path/to/blueprintFile.bed
		This option provides the path to a bed file which contains a list of contiguous genomic
		loci indicating the sites of known open chromatin sites, tagged with the number of cell
		types with open chromatin at that site. The format is as follows:

			chr1	1	10099	C#1	0
			chr1	10100	10330	#1	37
			chr1	10331	10344	C#2	0
			chr1	10345	10590	#2	4
			chr1	10591	16099	C#3	0

		where the columns respectively indicate the chromosome name, stie start position, site
		end position, element name, and number of cell types with open chromatin at that site.
		When this file is provided, PeaKDEck calculates signal-to-noise ratio, and calculates 
		the background probability distribution from sites selected from loci with no known 
		open chromatin.

	-npBack positiveInteger
		Sets the number of sites to randomly select to calculate the background probability 
		distribution. By default this is set to 50000 sites. -npBack must be a positive whole
		number.

	-sig probabilityValue
		Specifies the positive limit of the probability distrubtion for selecting the corrected
		read density for peak threshold. By default, -sig is set to 0.001. -sig must be a 
		positive number between 0 and 1.

	-PVAL ON|OFF
		Peaks are scored with the maximum corrected read density recorded during that peak by 
		default. setting -PVAL to ON converts this corrected read density to a probability value
		from the background probability distribution used to calculate the threshold. This
		value represents the probability that a corrected read density of that magnitude 
		belongs to the background probability distribution.

=head2 Top peak selection (-T I<peaks.bed>)
	
	Sorts peak bed files in descending order by corrected read density score. By default,
	the sorted peaks are printed in bed format to STDOUT. The target file must by in 
	bed format.

=head3 Mandatory settings

=head3 Optional settings
	
	-n positiveInteger
		This specifies the number of peaks to include in the resulting bed file, from the
		highest scoring peak downwards. By default, -n is set to ALL, and all the peaks 
		are printed to the output file. -n must be either 'ALL' or a positive whole number. 

=head1 DESCRIPTION

	PeaKDEck is a utility written in perl, mainly intended for use in the identification
	of peaks in mapped DNaseI-seq data. It also includes a set of utilities for 
	processing and manipulation of this data from the mapping stage forwards. It works
	on data in sam format.

	PeaKDEck selects a threshold read density for peak calling by constructing a 
	probability distribution of background read density scores using kernal density
	estimation. It selects a threshold by selecting a read density that is 'significantly' 
	outside this background probability distribution. All measurements of read density are 
	corrected for local background variation in signal intensity.

	PeaKDEck is also available as a standalone GUI application for all major platforms. 
	The GUI wrapper is written in Perl/Tk, and is available at 
	www.ccmp.ox.ac.uk/PeaKDEck.

=head1 FAQs

=head2 What are the system requirements?

	The command line and GUI PeaKDEck applications have been tested on OSX (Mountain Lion),
	Ubuntu 12.04 LTS, Windows XP and Windows 7. The system requirements are largely dependent on 
	the size of data files being used. We recommend at least 4GB memory for basic use with 
	small data files. For the numerical sorting of sam files, ~(file size * 2.5) free memory 
	is required for efficient sorting. For the command line applications, we recommend Perl 
	v5.12 or later. On Windows, PeaKDEck was tested with Strawberry Perl. 

=head2 How do I install PeaKDEck?

	PeaKDEck GUI: On Windows, PeaKDEck should run without the need to install Perl, or
	orther additional softwares. On Linux and OSX platforms, the X Window System (X11 or 
	XQuartz) must be installed.

	PeaKDEck command line: to use the PeaKDEck command line application, Perl must be 
	installed on your computer. We recommend Perl v5.12 or later. On Linux and OSX platforms,
	no other software is required to run the command line application. On the Windows 
	platform, a pseudorandom number generating module (Math::Random::MT) is required, and 
	is available through CPAN for Strawbery Perl users, and PPM (Math-Random-MT) for  
	ActiveState users. 

=head2 Which short read file formats does PeaKDEck work with?

	At present, PeaKDEck only works with files in the SAM format (see 
	samtools.sourceforge.net/SAMv1.pdf‎ for details). 

=head2 Which application should PeaKDEck be opened with?

	On the OSX platform, after launching the PeaKDEck GUI application (having installed 
	XQuartz), you may be prompted to choose an application with which to open PeaKDEck.
	PeaKDEck should be opened with the Terminal application (located at 
	/Applications/Utilities/Terminal.app).

=head2 Why is the GUI freezing?

	As yet, the PeaKDEck GUI is not a multi-threaded program. As such during data processing, 
	the GUI may appear frozen or unresponsive. Particularly on the Windows platform. For now,
	this is expected behaviour. The GUI will refresh when new status updates are available,
	and will return to full responsiveness when data processing has finished.

=head1 EXAMPLES

	see SYNOPSIS above.

=head1 CAVEATS

=head1 AUTHOR
	
	michael mccarthy
	michaeltmccarthy@gmail.com

=head1 ACKNOWLEDGEMENTS

=head1 SEE ALSO

=cut

