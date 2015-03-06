#!/usr/bin/perl
use warnings;
use strict;

print "---------------------------------------------------------------------------------------------------------------------------------------------\n";
print "Author: Steven Robbins\n";
print "This version of the script takes a normalised or non-normalised, tab separated, qiime style  OTU table as a positional argument and returns\n";
print "an L6 style table, converting tag counts to percentages and removing taxa that were not above  0.5% abundance in any one sample.  It also\n";
print "returns tax strings which have been truncated to the lowest identified taxonomic level.\n";
print "---------------------------------------------------------------------------------------------------------------------------------------------\n";
print "\n";

my $infile = $ARGV[0]; #normalised OTU table
my $outfile = "uncolapsed_"."$ARGV[0]";

open (IN, "<$infile") or die "could not open infile: $infile";
open (OUT, ">$outfile") or die "could not open outfile: $outfile";

#array of read counts for each column
my @read_counts;
#keeps track of row #. Don't need first two rows.
my $row_counter = 0;
#my $c2_counter = 0;
my $line_counter = 0;
my %OTUs;
my %OTU_count;

#could add one loop that checks and outputs the number of columns

#get read count of sample/column 1
while (my $line = <IN>)
{
        chomp $line;

        my @array = split(/\t/, $line);
        $row_counter++;
        my $column_counter = 0; #keeps track of the column when adding current cell value to read count

        #if at 3rd row or lower, add += each element to its appropriate @read_count element so that the array contains an ordered list of read counts
        if ($row_counter >= 3)
        {
                shift @array;
                pop @array;

                foreach my $element (@array)
                {
                        $read_counts[$column_counter] += $element;
                        $column_counter++;
                }
        }
}
close IN;

open (IN, "<$infile") or die "could not open infile: $infile";

while (my $line = <IN>)
{
        my @perc_array; #array containing perc abundance of tags
        chomp $line;
        $line =~ s/ //g;

        my @array = split(/\t/, $line);
        $line_counter++;

        #strip first line
        if ($line_counter == 1)
        {}

        #format and print header line for L6 file
        if ($line_counter == 2)
        {
                shift @array;
                unshift (@array, 'Taxon');
                pop @array;

                print OUT  join("\t", @array), "\n";
        }

        #process all data entries by converting to percentages and appending taxon string to beginning of line
        if ($line_counter >= 3)
        {
                #remove first column (OTU #) and last column (tax string) to add as key to %OTUs
                my $otu_id = shift @array;
                my $tax_string = pop @array;
                $tax_string =~ s/;s__$//;
                $tax_string =~ s/;g__$//;
                $tax_string =~ s/;f__$//;
                $tax_string =~ s/;o__$//;
                $tax_string =~ s/;c__$//;
                $tax_string =~ s/;p__$//;
                my $read_counts_position = 0; #will keep track of where I am in the read_counts array that countains the read counts for each sample

                #add data columns to perc_array as percentages
                foreach my $element (@array)
                {
                        my $perc_element = $element / $read_counts[$read_counts_position];
                        push (@perc_array, $perc_element);
                        $read_counts_position++;
                }

                #convert the array of data values into a string to be stored in the hash %OTUs
                my $array_string = join("\t", @perc_array);
#determine if the OTU is present in any sample at over 0.5%
                my $flag = 0;
                foreach my $element (@perc_array)
                {
                        #CHANGE THIS NUMBER IF YOU WANT TO CHANGE THE THRESHOLD FOR GETTING RID OF TAXA LESS THAN X% ABUNDANCE 
                        if($element > 0.005)
                        {
                                $flag=1;
                                        last;
                        }
                }

                if($flag==1)
                {

                        if (exists $OTUs{$tax_string})
                        {
                                $OTU_count{$tax_string}++;
                                $tax_string = "$tax_string"."-"."$OTU_count{$tax_string}";
                                $OTUs{$tax_string} = "$array_string"."  $otu_id";
                        }
                        else
                        {
                                $OTUs{$tax_string} = "$array_string"."  $otu_id";
                                $OTU_count{$tax_string} = 1;
                        }
                }
        }
}

my @hash_array = sort keys %OTUs;
foreach my $element (@hash_array)
{
        print OUT "$element\t$OTUs{$element}\n";
}

close IN;
close OUT;
exit;
