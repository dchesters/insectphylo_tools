
#
#	check_terminals.pl, used in pipeline for processing machine readable phylogenies
#		  
#    	Copyright (C) 2022-2024 Douglas Chesters
#
#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# 
# 
################################################################################################################ 
# 
# 
# 2022-06-27: 	Script start. 
#		after a bunch of different scripts each for processing a specific published newick tree,
# 		this script checks the result and outputs the terminals which conform to various specifications,
# 		output list is used for pruning, and making the standardized version of the tree.
# 2023-05-19:	Several 2letter species names in butterflies, now allowed for.
# 2023-08-22:	lists terminals that look like subfams / fams
# 2024-09-09:	Added licence for upload to github.
# 
# 
# 
# 
# 
# 
# 
##########################################################################################

$file = $ARGV[0];  # input
$file2 = $ARGV[1]; # outputs

open(IN, $file) || die "\nerror cant open file ($file)\n";
while(my $line = <IN>)
	{
	if($line =~ /\(/)
		{
		$tree = $line;

		
		};
	};
close IN;

open(OUT, ">$file2") || die "";
open(OUT2, ">$file2.tre") || die "";
open(OUT3, ">$file2.prune_from_file") || die "";

$treeCOPY = $tree;
$treeCOPY2 = $tree;


while($tree =~ s/([\(\,])([A-Z][^\:\,\)]+)([\:\,\)])/$1$3/)
	{
	my $terminal = $2;

	my $genus = "NA";
	if($terminal =~ /^([A-Z][a-z]+)/)
		{
		$genus = $1;
		if($genus =~ /^[A-Z][a-z]+i[dn]ae$/){$families{$terminal}=1};
		}else{print "\nWARNING, couldnt parse genus.\n"};
	$genus_counts{$genus}++;

	if($check_repeat{$terminal} == 1)
		{
		$counter++;
		print "\nWARNING, repeat terminal:$terminal\n";$count_duplicate_terminals++;
		$treeCOPY =~ s/([\(\,])($terminal)([\:\,\)])/$1 . $terminal . $counter . $3/e; $store_duplicate_terminal_IDs{ $terminal . $counter } =1;

		};
	$check_repeat{$terminal} = 1;

	if($terminal =~ /^[A-Z][a-z]+_[a-z]{3,30}$/)
		{
		print OUT ">$terminal\n"; # expected binomial
		}else{

		if($terminal =~ /_[a-z][a-z]$/)
			{
			# 2 lowercase letters in species field

			if($terminal =~ /_sp/)
				{
				# species ambiguous: sp.

				$non_binomial_terminals{$terminal} = $genus; 
				$non_binomial_count++; print "terminal:$terminal\tNON BINOMIAL\n";
				}else{
				# 2 letter but not sp., must be a species name of 2 letters
				print OUT ">$terminal\n";
				};

			}else{
			# not an expected binomial, not 2 letter, must be something weird like a code
			$non_binomial_terminals{$terminal} = $genus; 
			$non_binomial_count++; print "terminal:$terminal\tNON BINOMIAL\n";
			};

		};
#	$binomials_current++;$all_binomials{$binomial} = 1;

	};

print "
non_binomial_count:$non_binomial_count
";
my @non_binoms = keys %non_binomial_terminals; @non_binoms = sort @non_binoms;
print "\n\n";
foreach my $terminal(@non_binoms)
	{
	$belongs_to_genus 	= $non_binomial_terminals{$terminal};
	my $count 		= $genus_counts{$belongs_to_genus};

	if($count >= 2)
		{
		print "\t$terminal belongs to $belongs_to_genus which has congeners ($count), so not needed\n";
		$store_prune_IDs{$terminal} = 1;
		}else{
		$count_genus_unique_species_ambig++;print "although terminal $terminal is species ambiguous, it is genus unique so can be used\n";
		print OUT ">$terminal\n";
		};
	};
close OUT;

print OUT2 "$treeCOPY\n";
close OUT2;

$tree =~ s/\:[\d\.]+//g;




print "
count_genus_unique_species_ambig:$count_genus_unique_species_ambig
remaining string: $tree
";


my @duplicated_terminal_IDs = keys %store_duplicate_terminal_IDs;@duplicated_terminal_IDs = sort @duplicated_terminal_IDs;
if($count_duplicate_terminals >= 1)
	{
	open(DUPLICATED, ">Duplicated_Terminal_IDs") || die "";
	my $duplicated_list;
	foreach my $ID(@duplicated_terminal_IDs)
		{
		if($ID =~ /[\w\d]/){print DUPLICATED "$ID\n";$duplicated_list .= "$ID "};
		};
	close DUPLICATED;
	print "\nwarning, there are $count_duplicate_terminals duplicated terminal IDs in your tree. \n";
	print "have appended numbers to duplcates and written to new tree file:$file2.tre\n";
	print "list of these is written to file Duplicated_Terminal_IDs\n";
	print "OR, here is a list that can be used in switch for prune_tree\n-start_prune_list $duplicated_list-end_prune_list\n\n";

	};



while($treeCOPY2 =~ s/(\([^\,\(\)]+\,[^\,\(\)]+\))/internalnode/)
	{
	my $matched = $1; # print "matched:$matched\n";
	};

#print "\ntreeCOPY2:$treeCOPY2\n";

if(length($treeCOPY2) <= 20)
	{
	print "\ntree is balanced.\n"
	}else{
	print "\n\nWARNING, tree appears to be unbalanced.\n\n";
	};





my @pruneIDs = keys %store_prune_IDs;@pruneIDs = sort @pruneIDs;
foreach my $term(@pruneIDs)
	{
	print "uninformative:$term\n";
	print OUT3 "$term\n";
	};

my @fams = keys %families;@fams = sort @fams;
print "following terminals look like family/subfam level names:\n";
foreach my $term(@fams)
	{
	print "term:$term\n";
	if($store_prune_IDs{$term} == 1)
		{
		}else{
		print OUT3 "$term\n";
		};
	};
close OUT3;


