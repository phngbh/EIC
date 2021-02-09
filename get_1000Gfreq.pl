#One needs to install Ensembl API before using this script

use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

# Load the registry from db
$registry->load_registry_from_db(
  -host => 'ensembldb.ensembl.org',
  -user => 'anonymous'
);


my $va_adaptor = $registry->get_adaptor('human','variation','variation');
$va_adaptor->db->use_vcf(1);


my $filename = '100snps.txt';
open(my $fh, '<:encoding(UTF-8)', $filename)
  or die "Could not open file '$filename' $!";
 
my @snps = <$fh>;
chomp @snps;

foreach my $snps (@snps) {
  print "$snps\n";
}

foreach my $snps (@snps) {
  my $va = $va_adaptor->fetch_by_name($snps); #get the Variation from the database using the name
  my $genotypes = $va->get_all_PopulationGenotypes();
  foreach my $genotype (@{$genotypes}) {
   my $pop_name = $genotype->population->name;
   my $gt_name = $genotype->genotype_string;
   my $freq = $genotype->frequency;
    open (FH, '>>', 'genotypes2.txt') or die "Could not open file $!";
  print FH "$snps $pop_name $gt_name $freq\n";
  close (FH);
}
}

print('Write to genotype successfully!')

#foreach my $snps (@snps) {
#  my $vaa = $va_adaptor->fetch_by_name($snps); #get the Variation from the database using the name
#  my $alleles = $vaa->get_all_Alleles();

# foreach my $allele (@{$alleles}) {
#  next unless (defined $allele->population);
#  my $allele_string   = $allele->allele;
#  my $frequency       = $allele->frequency || 'NA';
#  my $population_name = $allele->population->name;
#  open (FH, '>>', '/home/ensembl/src/alleles2.txt') or die "Could not open file $!";
#  print FH "$snps $population_name $allele_string $frequency\n";
#  close (FH);
#}

#}

#print('Write to allele successfully!')


