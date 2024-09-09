

#
#	process_catalogue_of_life.pl, for parsing species from the CoL database
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
#################################################################################


$in = $ARGV[0];
$out = $ARGV[1];
if($in eq $out){die "\ncommand error.\n"};

$job=2;

if($job==2)
	{
	COL_to_taxontable();exit;
	};


open(IN, $in) || die "\ncant open in $in\n";
open(OUT, ">$out") || die "\ncant open OUT $out\n";
print "opened $in ...\n";
$in_insects=0;
while(my $line = <IN>)
	{
	$line =~ s/\n//;$line =~ s/\r//; $lines++;
	if($line =~ /\[class\]/ && $in_insects == 1){$in_insects++};
	if($line =~ /Insecta\s+\[class\]/){$in_insects++};
	if($in_insects == 1)
		{
		$insect_lines++;
		if($line =~ /^\s+\*/)
			{
			$synonyms_rm++
			}else{
			print OUT "$line\n"
			};
		};
	};
close IN;
close OUT;


print "
all lines:$lines
insect_lines:$insect_lines
synonyms_rm:$synonyms_rm

";



#####################################################################################################

sub COL_to_taxontable
{

open(IN, $in) || die "\ncant open in $in\n";
open(OUT, ">$out") || die "\ncant open OUT $out\n";
print "opened $in ...\n";
$in_insects=0;
while(my $line = <IN>)
	{
	$line =~ s/\n//;$line =~ s/\r//; $lines++; # print "$line\n";
	if($line =~ /^(\s+)(\S.+)\s\[(\w+)\]/)
		{
		my $step_string = $1; my $taxon_string = $2; my $rank = $3; $ranks{$rank}++;
		my $step_length = length($step_string); # print "step_length:$step_length\n";

		my $binomial;my $skip_line=0;

	#	if($rank eq "epifamily" || $rank eq "family" || $rank eq "genus" || $rank eq "infraorder" || $rank eq "order" || 
	#	$rank eq "parvorder" || $rank eq "species" || $rank eq "subfamily" || $rank eq "subgenus" || $rank eq "class" || 
	#	$rank eq "suborder" || $rank eq "subtribe" ||  $rank eq "superfamily" || $rank eq "supertribe" || $rank eq "tribe")		
		unless($rank eq "subspecies" || $rank eq "infraspecific_name")
			{
			if($rank eq "species")
				{
				# parsing complex
				
				$taxon_string =~ s/[\-\[\]]//g;

				# Panscopus pallidus Buchanan, 1927
				if($taxon_string =~ /^([A-Z][a-z]+\s[a-z]+)\s[A-Z][a-z].+$/)
					{
					$binomial = $1;$taxon_string=$binomial;
					}elsif($taxon_string =~ /^([A-Z][a-z]+)\s\([A-Z][a-z]+\)\s([a-z]+)\s[A-Z][a-z].+$/)
					{	# Tanyproctus (Tanyproctus) sichuanicus Keith, 2007
					my $genusname = $1; my $speciesname = $2;$binomial=$genusname . "_" . $speciesname;$taxon_string=$binomial;
					}elsif($taxon_string =~ /^([A-Z][a-z]+\s[a-z]+)\s\(\S+\s\d{4}\)$/)
					{
					# :Bodilus zarudnyi (Frolov, 2001)
					$binomial = $1;$taxon_string=$binomial;
					}elsif($taxon_string =~ /^([A-Z][a-z]+\s[a-z]+)\s.+/)
					{
					$binomial = $1;$taxon_string=$binomial;
					}elsif($taxon_string =~ /^([A-Z][a-z]+)\s\([A-Z][a-z]+\)\s([a-z]+)\s.+/)
					{
					my $genusname = $1; my $speciesname = $2;$binomial=$genusname . "_" . $speciesname;$taxon_string=$binomial;
					}elsif($taxon_string =~ /^([A-Z][a-z]+\s[a-z]+)$/)
					{
					$binomial = $1;$taxon_string=$binomial;
					}elsif($taxon_string =~ /^([A-Z][a-z]+)\s\([A-Z][a-z]+\)\s([a-z]+)$/)
					{
					my $genusname = $1; my $speciesname = $2;$binomial=$genusname . "_" . $speciesname;$taxon_string=$binomial;
					}else{
					print "rank species:$taxon_string\n";$skip_line=1;
					# unconventially named stuff, and errors, and extict, ignore these.
					};
				}else{
				# not species, parsing 'easier'
				if($taxon_string =~ /^[A-Z][a-z]+$/)
					{
					# print "rank non-species, expected format:$taxon_string\n";
					# name only
					}else{
					# print "rank non-species, complex:$taxon_string\n";
					# taxon name plus persons name plus year
					if($taxon_string =~ s/^([A-Z][a-z]+)\s\S.+\s\d{4}$/$1/)
						{
						}elsif($taxon_string =~ s/^([A-Z][a-z]+)\s\S.+$/$1/)
						{
						
						}else{
						# print "!rank non-species:$taxon_string\n";
						# discard these and decendents; 
						# mostly exinct stuff, plus unconventionally named, plus errors.
						# pecies:Catopsylla†
						$skip_line = 1;
						};
					};
				};			

			$taxon_string =~ s/\s+/_/g;
			if($rank eq "species")
				{
				if($all_species{$taxon_string} == 1)
					{$skip_line = 1;print "warning, duplicate sp:$taxon_string\n"};
				$all_species{$taxon_string} = 1;
				};
			if($rank eq "genus")
				{
				if($all_genera{$taxon_string} == 1)
					{
				#	$skip_line = 1;
					print "\twarning, duplicate genus:$taxon_string\n"
					};
				$all_genera{$taxon_string} = 1;
				};

			if($skip_line == 0)
				{
				$current_lineage{$step_length} = $lines; $code_is_taxon{$lines} = $taxon_string;
				my $parent_tax = $current_lineage{ ( $step_length -2) };
				my $parent_name = $code_is_taxon{$parent_tax};
				# print "parent of $taxon_string is $parent_tax\n";
				if($lines == 1){$parent_tax = 0;$parent_name = "Root"};
				print OUT "$lines\t$parent_tax\t$taxon_string\t$parent_name\t$rank\n";
				};

			}; # only main ranks



		}else{
		print "not parsed:$line\n";
		};

	};

close IN;
close OUT;

print "
all lines:$lines
";	
	
my @ranks_array = keys %ranks;@ranks_array = sort @ranks_array;
print "@ranks_array\n";
# class epifamily family form genus infraorder infraspecific_name nanorder order other parvorder 
# species subfamily subgenus suborder subspecies subtribe superfamily supertribe tribe variety


};

#####################################################################################################













