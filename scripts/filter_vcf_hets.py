#!/usr/bin/env python

# This python script will convert a VCF file to hapmap files suitable for use in TASSEL
## or rrBLUP. The TASSEL version encodes genotypes by their nucleotides (i.e. AA, AG, or GG), 
## while the rrBLUP version encodes homozygotes for the reference allele as 1, heterozygotes at 0, and
## homozygous alternates as -1

# Import modules
import argparse # To get the arguments
import sys



#####
# Define the arguments
#####

# Description
DESC = """A python script to filter heterozygous genotypes in a VCF file.\n\n"""

# Argument parser
parser = argparse.ArgumentParser(description=DESC, add_help=True)

# Add arguments
# Input VCF file
parser.add_argument('-i',
  '--vcf',
  metavar = 'VCFIN',
  help = 'Input VCF file',
  required = True)
# Output file name
parser.add_argument('-o',
  '--outfile',
  metavar = 'OUTFILE',
  help = 'Output file basename (i.e. no extension)',
  required = True)
  
# Minimum allele depth of het to keep genotype call
parser.add_argument('-d',
  '--mindepth',
  metavar = 'DEPTH',
  help = 'Minimum depth of both alleles to retain a heterozygous genotype call.',
  required = True)
  
# Parse the arguments
args = parser.parse_args()

# Convert mindepth to int
mind = int(args.mindepth)


#####
# Execute program functions
#####


with open(args.vcf, 'r') as VCF:
  
  # Create the output filename
  filename = str(args.outfile) + '.vcf'
  # Open handle for writing file
  handle = open(filename, 'w')
  
  # Start reading through the vcf file
  for line in VCF:
    # Print all metadata lines
    if line.startswith('#'):
      handle.write(line)

      
    # Now we have the metadata and sample  names; let's move on to the genotype data
    else:
      # Split the line by tabs
      tmp = line.strip().split('\t')
      
      # Everything except the genotype calls
      meta = tmp[:9]
  
      # The genotype data
      genotypes = tmp[9:]
      # Empty list to store genotypes
      geno_out = []
  
      # Iterate over genotypes
      for g in genotypes:
        # The genotype string is separated by ':'
        g1 = g.split(":")
        # The first element of the genotype string is the genotype call
        call = g1[0]
        
        # If the call is heterozygous, get the depth
        if call == "0/1" or call == "1/0":
          
          # Allele depth and split by comma
          ad = g1[1].split(",")
          
          # Test if depth is greater than the provided minimum
          test = [int(d) < mind for d in ad]
          
          # If not all are true, make the genotype call missing
          if (any(test)):
            g1[0] = "./."
          
          # Collapse g1
          g2 = ":".join(g1)
          
        else:
          g2 = g
          
        # Append the genotype list
        geno_out.append(g2)
        
      # Combine meta with geno_out
      tmp1 = meta + geno_out
      # Write this line
      handle.write('\t'.join(tmp1) + '\n')
      
  print("File was written as " + filename)
  # Close the handle
  handle.close()
