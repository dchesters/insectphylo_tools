
#
#	taxon_table_to_newick.pl, converts taxonomy table to newick format
#		  
#    	Copyright (C) 2023-2024 Douglas Chesters
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
#################################################################################
# 
# 
# 
# 2023-08-26:	Script start. 
#		Converts taxonomic heirachy in tablur format to newick.
# 		Most code is from relational_constraints.pl
# 2024-09-09:	Added license.
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
#################################################################################


$in		= $ARGV[0];
$out		= $ARGV[1];
$starting_node	= $ARGV[2];

if($in eq $out){die "\ncommand error\n"};

print "
SETTINGS
  in:$in
  out:$out
  starting_node:$starting_node
";


# sub store_tax_heirarchy stores the following items:
# 	ncbi_nodes{}{rank_code}; $ncbi_nodes{}{rank}; $ncbi_nodes{}{parent}; 
#	$ncbi_nodes{}{child_nodes}; $ncbi_nodes{}{name}

#######################
store_tax_heirarchy();#
#######################


#################################################

# traverse nodes starts at first parent node in the input file, which should be the root.
# the input file was made with a root node user specified.


# sub traverse_nodes uses:
#  ncbi_nodes{}{child_nodes}; $ncbi_nodes{}{rank};$ncbi_nodes{}{name};$ncbi_nodes{}{rank_code};
# and writes:
#  complete_lineage_for_this_species ; ncbi_taxnumber_for_taxname	

open(LOG, ">taxon_table_to_newick_LOG");

$starting_name = $ncbi_nodes{$starting_node}{name};$starting_name_ff = $starting_name;$starting_name_ff =~ s/\s+/_/g;
print "name for the starting node ($starting_node)\nwhich has been specified ($starting_name)\n";
$starting_name =~ s/^(\w).+$/$1/;

$as_newick = "($starting_node)";

#################################################
traverse_nodes($starting_node , $starting_name);#
#################################################


my @ranks_array = keys %how_many;
print "heirarchy traversed:
	nodes_traversed:$nodes_traversed
	ranks encountered:$#ranks_array
	duplicate_species_removed:$duplicate_species_removed

";
close LOG;

$unconverted=0;
open(OUT2, ">$out.taxlist") || die "\nerror 79\n";
open(OUT3, ">$out.lineages") || die "\nerror 82\n";

my @all_terminals2 = keys %all_terminals;@all_terminals2 = sort @all_terminals2;
foreach my $terminal(@all_terminals2) # $terminal is just a number/ID, need to convert to name (includes species and higher taxa)
	{
	$ounter++;
	my $convert_to = "NA_";my $convert_to2 = "NA_";
	if($ncbi_nodes{$terminal}{name} =~ /\w/)
		{
		$convert_to = $ncbi_nodes{$terminal}{name};$convert_to2 = $ncbi_nodes{$terminal}{name};
		$converted++;
		if($terminal_encountered_already{$convert_to} == 1)
			{
			$duplicates_found++;$convert_to2 .= $duplicates_found;
			if($duplicates_found < 100)
				{
				print "warning, duplicate written to newick:$convert_to\n"
				}elsif($duplicates_found == 100)
				{
				print "wont bother printing more error messgs\n\n";
				}
			};
		$terminal_encountered_already{$convert_to} = 1;
		}else{
		$unconverted++;$convert_to .= $unconverted;
		};	
	if($as_newick =~ s/([\(\,])$terminal([\)\,])/$1$convert_to2$2/)
		{
		}else{
		$swap_errors++; # print "\twarning, didnt switch ID ($terminal) to name ($convert_to2) in newick\n";
		};

	my $lineage = "NA";
	if($complete_lineage_for_this_species{$convert_to2} =~ /\w/){ $lineage = $complete_lineage_for_this_species{$convert_to2}  };
	print OUT3 "$convert_to2\t$lineage\n";


	if($ounter =~ /0000$/)
		{
		print "$converted, code:$terminal, convert_to:$convert_to\n";
		}
	
	};

OUT3;

print "
duplicates_found:$duplicates_found out of:$converted
swap_errors:$swap_errors
";
$x=1;
if($x== 2)
{
while($as_newick =~ /([\(\,])(\d+)([\)\,])/)
	{
	my $char1 = $1;my $code = $2; my $char2 = $3;
	my $convert_to = "NA_";
	if($ncbi_nodes{$code}{name} =~ /\w/)
		{	
		$convert_to = $ncbi_nodes{$code}{name};
		$converted++;
		}else{
		$unconverted++;$convert_to .= $unconverted;
		};	
	$as_newick =~ s/([\(\,])(\d+)([\)\,])/$char1$convert_to$char2/;

	if($converted =~ /00$/)
		{
		print "$converted, code:$code, convert_to:$convert_to\n";
		}
	print OUT2 "$convert_to\n";
	};

};


print "
converted:$converted
unconverted:$unconverted
";
open(OUT, ">$out") || die "\nerror 49\n";
print OUT "$as_newick;\n";
close OUT;
close OUT2;

print "

Newick written to out file name $out
FIN.

";





#####################################################################################################
#
#
#
#####################################################################################################


sub store_tax_heirarchy
{
print "
reading users taxonomic file $taxon_table
";

$taxon_table = $in;
open(TAXTABLE, $taxon_table) || die "\ncant find taxon table ($taxon_table).\n";
my $tax_table_line_count=0;
while (my $line = <TAXTABLE>)
	{

# 	"$child\t" . 		# child node ID
#	"$current_node\t" . 	# parent node id
#	"$name_string\t" . 	# child tax
#	"$parentname\t" . 	# parent tax
#	"$rank\n"; 		# child rank

	$line =~ s/\n//;$line =~ s/\r//;
	if($line =~ /^(.+)\t(.+)\t(.+)\t(.+)\t(.+)/)
		{
		my $child_ID = $1;my $parent_ID = $2;my $child_tax = $3;my $parent_tax = $4;my $child_rank = $5;
		#print "child_ID:$child_ID parent_ID:$parent_ID child_tax:$child_tax parent_tax:$parent_tax child_rank:$child_rank\n";
		# child_ID:92482 parent_ID:92069 child_tax:Muscidora parent_tax:Trachyderina child_rank:genus
		# child_ID:92483 parent_ID:92482 child_tax:Muscidora_alutacea parent_tax:Muscidora child_rank:species
		# child_ID:92484 parent_ID:92482 child_tax:Muscidora_bezarki parent_tax:Muscidora child_rank:species


		$rank_hash{$child_rank}++;
		if($tax_table_line_count == 0)
			{
			$root_taxon_name = $parent_tax;
			# $starting_node = $parent_ID; 			# root node should be first parent node in the file
			$ncbi_nodes{$parent_ID}{name} = $parent_tax; 	# as names are assigned only for child nodes below, 
			}; 						# the root node needs name assigning by parent


#		my $tax_id = $1;my $parent_tax_id = $2;my $rank = $3;
#		print "tax_id:$tax_id parent_tax_id:$parent_tax_id rank:$rank\n";

		my $current_rankcode = 0;
		foreach my $rank_test("order","suborder","infraorder","series","superfamily","family","subfamily","tribe","subtribe","genus",
			"subgenus","species group","species","subspecies")
				{$current_rankcode++;
				if($child_rank eq $rank_test){$ncbi_nodes{$child_ID}{rank_code} = $current_rankcode}
				};

			$ncbi_nodes{$child_ID}{rank} = $child_rank;
			$ncbi_nodes{$child_ID}{parent} = $parent_ID;

			# 20231121 	as names are reused across ranks, eg subgenus name often same as its genus,		
			#		should include rank in node ID
			#		in fact better to use unique IDs and retrieve later name of each
			#		edit here using codes not names..

			if($ncbi_nodes{$parent_ID}{child_nodes} =~ /\t$child_ID\t/)
				{
				print "DUPLICATE 229:$child_ID,$child_tax\n";
				}else{
				$ncbi_nodes{$parent_ID}{child_nodes} .= "\t" . $child_ID . "\t";
				};

	
		$child_tax =~ s/[\(\)\,\[\]\'\#\&\/\:\.\-]/ /g;
		$child_tax =~ s/\s\s+/ /g;$child_tax =~ s/\s+$//;$child_tax =~ s/^\s+//;

		$child_tax =~ s/^Candidatus\s(\w+)$/$1/;

		$ncbi_nodes{$child_ID}{name} = $child_tax;


		$tax_table_line_count++;
		}else{
		print "cant parse line:$line\n";
		};

	};

close TAXTABLE;


print "taxon table has been read. 
	root taxon name:$root_taxon_name
 	node stored:$tax_table_line_count
";



}; # store_tax_heirarchy


#####################################################################################################
#
#
#
#####################################################################################################



#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################





sub traverse_nodes
{
my $current_node = $_[0];my $current_node_taxstring = $_[1];# $current_node_taxstring to be deprecated
$subcalls++;


 # for the current node, read the list of child nodes
my $child_nodes = $ncbi_nodes{$current_node}{child_nodes};$child_nodes =~ s/^\t//;$child_nodes =~ s/\t+$//; # NB these are CODES not NAMES
my @child_nodes_array = split(/\t+/, $child_nodes);
my $child_string_duplicate_species_removed = "";
foreach my $node7(@child_nodes_array)
	{
	if($node7 =~ /\d/) # ? only terminals have numbers?
		{$all_terminals{$node7}=1};	

	# here seems a good place to remove duplicate species	
	if($ncbi_nodes{$node7}{name} =~ /\w/)
		{
		my $name_of_ID = $ncbi_nodes{$node7}{name}; # print "node7:$node7 name_of_ID:$name_of_ID\n";
		if($name_of_ID =~ /^[A-Z][a-z]+_[a-z]+$/)
			{
			if($duplicate_species{$name_of_ID} == 1)
				{
				$duplicate_species_removed++;print "\n\tDUPLICATE SPECIES:$name_of_ID\n";
				}else{
				$child_string_duplicate_species_removed .= "	$node7	";
				};
			$duplicate_species{$name_of_ID} = 1;
			}else{
			$child_string_duplicate_species_removed .= "	$node7	";
			};

		}else{
		$child_string_duplicate_species_removed .= "	$node7	";
		};
	};

$nodes_traversed++;


 # seems mostly discontinued:
if($current_node == $starting_node)
	{
	my $rank = $ncbi_nodes{$starting_node}{rank};$rank =~ s/\s/_/g;$how_many{$rank}++;
#	my $current_node_taxstring = substr $current_node_taxstring, 0,1;
#	my $originalname = $ncbi_nodes{$starting_node}{name};$originalname =~ s/[\s\t]/_/g;
#	my $child_complete_lineage = $ncbi_nodes{$starting_node}{complete_lineage} . "$rank:$originalname ";
	};

my $parent_name_string = $ncbi_nodes{$current_node}{name};$parent_name_string =~ s/\s+/_/g;
my $parent_rank = $ncbi_nodes{$current_node}{rank};$parent_rank =~ s/\s+/_/g;

# this might contain duplicate species names:
#my $child_nodes_string = $child_nodes;
my $child_nodes_string = $child_string_duplicate_species_removed;

$child_nodes_string =~ s/^\t+//;$child_nodes_string =~ s/\t+$//;
$child_nodes_string =~ s/\t+/,/g;

print LOG "$current_node\t$parent_name_string\t$child_nodes_string\n";
if($subcalls =~ /0000$/){print "sub traverse_nodes calls $subcalls\tnodeID $current_node\tname $parent_name_string\n"};


# HERE building the newick string from the current traverse, at the moment comprising nodeIDs not names.
if($child_nodes_string =~ /\,/)
	{
	$as_newick =~ s/([\(\,])$current_node([\)\,])/$1($child_nodes_string)$parent_rank.$parent_name_string$2/;
	}elsif($child_nodes_string =~ /\d/)
	{
	$as_newick =~ s/([\(\,])$current_node([\)\,])/$1$child_nodes_string$2/;
	# print "\nunexpected 190, $child_nodes_string\n"
	};

foreach my $child(@child_nodes_array)
	{
	unless($current_node == $child) # one of the child nodes of the root node (1), is also 1 (NCBI system)
	{
	my $rank = $ncbi_nodes{$child}{rank};$rank =~ s/\s/_/g;$how_many{$rank}++;
	my $name_string = $ncbi_nodes{$child}{name};$name_string =~ s/\s+/_/g; ####################### sep2013 .. fixed mar2015
	my $rankcode = $ncbi_nodes{$child}{rank_code};$rank_codes{$name_string} = $rankcode;


	unless($ncbi_nodes{$current_node}{complete_lineage} =~ /\w/)#22oct2015:otherwise where shared taxon is same as basal node, wont be assigned 
		{$ncbi_nodes{$current_node}{complete_lineage}="root:$starting_name_ff "};

	my $child_complete_lineage = $ncbi_nodes{$current_node}{complete_lineage} . "$rank:$name_string ";#prev assigned $name


	# NOV 2018: user decides which ranks are constrained	
	my $constrain_this_node = 0; # print "\n\n";
	foreach my $user_constrain_rank ( @constrain_ranks )
		{
		if($rank eq $user_constrain_rank)
			{$constrain_this_node = 1};
	#	print "rank:$rank user:$user_constrain_rank constrain:$constrain_this_node\n";
		};
		$ncbi_nodes{$child}{complete_lineage} = $child_complete_lineage;


	#print "storing complete lineage for taxon:$name_string\n";	
	$complete_lineage_for_this_species{$name_string}=$child_complete_lineage;
#	print "child_complete_lineage:$child_complete_lineage\n";

	$ncbi_tax_number_for_this_species{$name_string}=$child;
	my $originalname = $ncbi_nodes{$child}{name};$originalname =~ s/[\s\t]/_/g;
	$ncbi_taxnumber_for_taxname{$originalname} = $child;
	my $name_assignment_to_taxnumber = "";

	if($ignore_subspecies == 1)
		{
		if ( $rank eq "subspecies" || $rank eq "varietas"   || $nodes{$current_node}{complete_lineage} =~ / species:/)# 
			{
			my $parentoriginalname = $ncbi_nodes{$current_node}{name};
			if($parentoriginalname =~ /^[A-Z][a-z]+\s[a-z]+$/)
				{
				$parentoriginalname =~ s/[\s\t]/_/g;
				$name_assignment_to_taxnumber = $parentoriginalname;
			#	print "node:$child named:$originalname rank:$rank appears to be subspecfic\n";
			#	print "\tassigning parent name instead:$parentoriginalname\n\n";
				}else{$name_assignment_to_taxnumber = $originalname}
			}else{
			$name_assignment_to_taxnumber = $originalname
			}

		}else{
		$name_assignment_to_taxnumber = $originalname
		}

	#print "$name_assignment_to_taxnumber $child $rank $child_complete_lineage\n";


		###########################################
		traverse_nodes($child , $child_taxstring);#
		###########################################
	}}


	
}#sub traverse_nodes





#####################################################################################################
#
#
#
#####################################################################################################




