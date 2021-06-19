#! /usr/bin/perl -w

use FileHandle;

# Script to run batch jobs for the cteqX global fit project
# Author : Peter Monaghan - 19th November 2008 
# This script uses the most up to date syntax for submitting
# jobs to the Auger batch farm, including the project code
# for the cteqX group and the track name.

# Further information can be found on the JLab Scientific
# Computing webpages; look for a link to the
# "Auger_XML_Configuration_File_Format_Guide".

# When executing this script, a single input file name
# is required on the command line. This input file will contain 
# executable and input file names for each job to be submitted.
# A separate xml file is created for each individual job, saved 
# and then automatically submitted to the batch farm with 
# the jsub command

my $rfile=$ARGV[0];   # input file name

# Setting some tags for use later
# CHANGE the directory paths as necessary for your system

my $dir="/group/cteqX/peter/fitpack/input"; 
my $runfile="$dir/$rfile";
print "$runfile\n";

# Create a generic file tag for the output files
my $outf=new FileHandle;

# Open the input file

open(RUN_FILE,$runfile) or die("Can't open file $rfile");
my $Nrun;  

# Using a while loop to read in the input file line by line.
# A separate xml file is created for each input line.

while(<RUN_FILE>){
    
    # split each line into its components
    # i.e. executable, input, accept or reject
    my @lines=split(/\#/,$_,3);
    my @line =split(/ /,$lines[0],3);
    $line[1]=~ s/\s//g;  # remove spaces 
    $line[2]=~ s/\s//g;  # remove spaces 
    $Nrun = $.;          # get line number


    my $jsubfile;
    if($line[2] eq "accept") {
	$jsubfile="jsub\_$line[0]\_$Nrun.xml"; #create xml filename
	print "$jsubfile\n";              # print out filename
	print "line number = $Nrun \n";
	print "executable = $line[0] \n";
	print "input file = $line[1] \n";
	print "status = $line[2] \n";

	# Create the output xml file
	$outf->open(">$jsubfile") or die("Can't open file $jsubfile for writing");

	# Write out the necessary syntax for the xml submission file
	print $outf "<Request>\n";
	
	# CHANGE email address as necessary
	print $outf "   <Email email=\"peter\@jlab.org\" request=\"false\" job=\"true\"/>\n";
	# our project name is simply cteqX
	print $outf "   <Project name=\"cteqX\"/>\n";
	# Track name is a batch farm parameter
	# 'reconstruction' is fine for our purposes
	print $outf "   <Track name=\"reconstruction\"/>\n";
	# Name is whatever you want - CHANGE as necessary
	print $outf "   <Name name=\"cteqfit\"/>\n";
	print $outf "\n";  # inserts a blank line

	# Now we are going to specify any directory paths 
	# which any of our codes or jobs might require

	# ADD or DELETE or CHANGE for your own particular case 
	# CHANGE the paths for your own set of folders 
	# where you are running the code from!!!

	print $outf "   <Variable name=\"dir\" value=\"file:/group/cteqX/peter/\"/>\n";
	print $outf "   <Variable name=\"fitdir\" value=\"file:/group/cteqX/peter/fitpack/fitting/\"/>\n";
	print $outf "   <Variable name=\"datadir\" value=\"file:/group/cteqX/peter/fitpack/data/\"/>\n";
	print $outf "   <Variable name=\"theorydir\" value=\"file:/group/cteqX/peter/fitpack/theory/\"/>\n";
	print $outf "   <Variable name=\"utildir\" value=\"file:/group/cteqX/peter/fitpack/util/\"/>\n";
	print $outf "   <Variable name=\"nucldir\" value=\"file:/group/cteqX/peter/fitpack/nucl/\"/>\n";
	print $outf "   <Variable name=\"outdir\" value=\"file:/group/cteqX/peter/output/\"/>\n";
	print $outf "\n";

	# Next we specify the source files that should be copied 
	# from your working directory to the local batch farm computer

	print $outf "   <Input src=\"\${fitdir}/$line[0]\" dest=\"$line[0]\"/>\n";
	print $outf "   <Input src=\"\${fitdir}/$line[1]\" dest=\"$line[1]\"/>\n";
	print $outf "\n";

	# Specify amount of memory required
	# Too high a number can cause jobs to sit pending on 
	# the batch farm waiting for a machine with enough 
	# memory to become available
	
	print $outf "   <Memory space=\"200\" unit=\"MB\"/>\n";

	# Now specify the actual command(s) to execute to run 
	# the fitting code with the specified input file.
	# This defines the actual 'job'.

	print $outf "   <Job>\n";

	# Now we explicitly set some paths so that the 
	# code will work correctly on the local batch farm computer
	# REMEMBER TO SET THE DIRECTORY PATH APPROPRIATE FOR YOU

	print $outf "   <Command><![CDATA[\n";
	print $outf "      setenv cteqx_dat /group/cteqX/peter/fitpack/data/\n";
	print $outf "      echo \$cteqx_dat\n";

	print $outf "      echo $line[1] | ./$line[0]\n";
	print $outf "   ]]></Command>\n";
	print $outf "\n";

	# Include the relocation of the output files
	# assuming the names are 'fixed'
	# Send the output to whatever you want in YOUR directory

	print $outf "   <Output src=\"*.out\" dest=\"\${outdir}/test/*.out\"/>\n";
	print $outf "   <Output src=\"*.pdf\" dest=\"\${outdir}/test/*.pdf\"/>\n";
	# Capture the standard output to a file in your directory
	print $outf "   <Stdout dest=\"\${outdir}/summary/output_$Nrun.log\"/>\n";
	# Capture the standard error to a file in your directory
	print $outf "   <Stderr dest=\"\${outdir}/summary/outerr_$Nrun.log\"/>\n";

	print $outf "\n";

	print $outf "   </Job>\n";
	print $outf "\n";
	print $outf "</Request>\n";

	# Close the created xml file
	$outf->close;
    
	# Create the command to submit this xml file to the batch farm
	my $sub_cmd="jsub -xml $jsubfile";

	system($sub_cmd);   # execute the submission command
	
	# save the created xml file 
	# move it wherever you want to save it
	system("mv $jsubfile /group/cteqX/peter/output/");  

   }  # end of this line; continue to the start of the loop for the next line 

} # end of the original while loop

print "end of input file\n";
