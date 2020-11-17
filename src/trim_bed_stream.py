#!/usr/bin/env python3


"""
Program: bed_trim.py
Purpose: remove from bed file rows with bad co-ordinates.
Notes:
  o  bbi-sciatac-analyze/src/get_unique_fragments.py identifies a fragments site and
     writes the site co-ordinates as site_loc +/- n bases, where n=1 by default. Bad
     entries can result.
  o  this program skips bad entries
"""


import sys


def read_chromosome_size( file ):
  chromosome_size_dict = {}
  for line in file:
    fields = line.split()
    if( chromosome_size_dict.get( fields[0] ) ):
      print( 'Error: duplicate chromosome name \'%s\' in chromosome size file.' % ( fields[0] ), file=sys.stderr )
      sys.exit( -1 )
    chromosome_size_dict.setdefault( fields[0], int( fields[1] ) )
  return( chromosome_size_dict )


def bed_trim( file_in, chromosome_size_dict, file_out ):
  chromosome_name = None
  for line in file_in:
    fields = line.split()
    if( fields[0] != chromosome_name ):
      chromosome_name = fields[0]
      chromosome_end  = chromosome_size_dict.get( fields[0] )
      if( not chromosome_end ):
        print( 'Error: unknown chromosome name \'%s\'' % ( fields[0] ), file=sys.stderr )
        sys.exit( -1 )
    if( int( fields[1] ) < 0 or int( fields[2] ) >= chromosome_end ):
      print( 'remove %s' % ( line ), file=sys.stderr )
      continue
    print( '%s' % ( '\t'.join( fields ) ) , file=file_out )
  return( 0 )


if __name__ == "__main__":


  chromosome_size_filename = sys.argv[1]
  chromosome_size_dict = read_chromosome_size( open( chromosome_size_filename, 'rt' ) )
  bed_trim( sys.stdin, chromosome_size_dict, sys.stdout )

