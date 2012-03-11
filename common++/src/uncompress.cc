/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkCommon.h>

#ifndef SACREDMEM
#define SACREDMEM	0
#endif

#ifndef USERMEM
# define USERMEM 	450000
#endif

#ifdef USERMEM
# if USERMEM >= (433484+SACREDMEM)
#  define PBITS	16
# else
#  if USERMEM >= (229600+SACREDMEM)
#   define PBITS	15
#  else
#   if USERMEM >= (127536+SACREDMEM)
#    define PBITS	14
#   else
#    if USERMEM >= (73464+SACREDMEM)
#     define PBITS	13
#    else
#     define PBITS	12
#    endif
#   endif
#  endif
# endif
# undef USERMEM
#endif

#ifdef PBITS
# ifndef BITS
#  define BITS PBITS
# endif
#endif

#if BITS == 16
# define HSIZE	69001
#endif
#if BITS == 15
# define HSIZE	35023
#endif
#if BITS == 14
# define HSIZE	18013
#endif
#if BITS == 13
# define HSIZE	9001
#endif
#if BITS <= 12
# define HSIZE	5003
#endif

#define FIRST	257
#define	CLEAR	256

#define BIT_MASK	0x1f
#define BLOCK_MASK	0x80

#define INIT_BITS       9

#if BITS > 15
typedef long int	code_int;
#else
typedef int		code_int;
#endif

typedef long int	count_int;
typedef	unsigned char	char_type;

int n_bits;				/* number of bits/code */
int maxbits = BITS;			/* user settable max # bits/code */
code_int maxcode;			/* maximum code, given n_bits */
code_int maxmaxcode = 1L << BITS;	/* should NEVER generate this code */

#ifdef COMPATIBLE		        /* But wrong! */
#define MAXCODE(n_bits)	                (1L << (n_bits) - 1)
#else
#define MAXCODE(n_bits)	                ((1L << (n_bits)) - 1)
#endif

count_int htab [HSIZE];
unsigned short codetab [HSIZE];
code_int hsize = HSIZE;			/* for dynamic table sizing */
count_int fsize;

#define tab_suffixof(i)	                ((char_type *)(htab))[i]
#define de_stack		        ((char_type *)&tab_suffixof(1<<BITS))

code_int free_ent = 0;			/* first unused entry */

int block_compress = BLOCK_MASK;
int clear_flg = 0;

char_type lmask[9] = {0xff, 0xfe, 0xfc, 0xf8, 0xf0, 0xe0, 0xc0, 0x80, 0x00};
char_type rmask[9] = {0x00, 0x01, 0x03, 0x07, 0x0f, 0x1f, 0x3f, 0x7f, 0xff};

#ifdef HAS_ZLIB

code_int getcode(gzFile file, int fflg)
{
  code_int code;
  static int offset = 0, size = 0;
  static unsigned char buf[BITS];
  int r_off, bits;
  unsigned char *bp = buf;

  if (fflg) {
    offset=0;size=0;
  }

  if ( clear_flg > 0 || offset >= size || free_ent > maxcode ) {
    if ( free_ent > maxcode ) {
      n_bits++;
      if ( n_bits == maxbits ) {
        maxcode = maxmaxcode;	/* won't get any bigger now */
      } else {
        maxcode = MAXCODE(n_bits);
      }
    }
    if ( clear_flg > 0) {
      maxcode = MAXCODE (n_bits = INIT_BITS);
      clear_flg = 0;
    }
    size = gzread(file, (char *)buf, n_bits);
    if ( size <= 0 ) {
      return -1;			/* end of file */
    }
    offset = 0;
    size = (size << 3) - (n_bits - 1);
  }

  r_off = offset;
  bits = n_bits;
  bp += (r_off >> 3);   /*------ Get to the first byte.------- */
  r_off &= 7;
  code = (*bp++ >> r_off);  /* Get first part (low order bits) */
  bits -= (8 - r_off);
  r_off = 8 - r_off;	     /*---- now, offset into code word ----*/

  /* Get any 8 bit parts in the middle (<=1 for up to 16 bits). */
  if ( bits >= 8 ) {
    code |= *bp++ << r_off;
    r_off += 8;
    bits -= 8;
  }

  /*------ high order bits. -------*/
  code |= (*bp & rmask[bits]) << r_off;
  offset += n_bits;

  return code;
}

#else

code_int getcode(FILE *file, int fflg)
{
  code_int code;
  static int offset = 0, size = 0;
  static unsigned char buf[BITS];
  int r_off, bits;
  unsigned char *bp = buf;

  if (fflg) {
    offset=0;size=0;
  }

  if ( clear_flg > 0 || offset >= size || free_ent > maxcode ) {
    if ( free_ent > maxcode ) {
      n_bits++;
      if ( n_bits == maxbits ) {
        maxcode = maxmaxcode;	/* won't get any bigger now */
      } else {
        maxcode = MAXCODE(n_bits);
      }
    }
    if ( clear_flg > 0) {
      maxcode = MAXCODE (n_bits = INIT_BITS);
      clear_flg = 0;
    }
    size = fread(buf, 1, n_bits, file);
    if ( size <= 0 ) {
      return -1;			/* end of file */
    }
    offset = 0;
    size = (size << 3) - (n_bits - 1);
  }

  r_off = offset;
  bits = n_bits;
  bp += (r_off >> 3);   /*------ Get to the first byte.------- */
  r_off &= 7;
  code = (*bp++ >> r_off);  /* Get first part (low order bits) */
  bits -= (8 - r_off);
  r_off = 8 - r_off;	     /*---- now, offset into code word ----*/

  /* Get any 8 bit parts in the middle (<=1 for up to 16 bits). */
  if ( bits >= 8 ) {
    code |= *bp++ << r_off;
    r_off += 8;
    bits -= 8;
  }

  /*------ high order bits. -------*/
  code |= (*bp & rmask[bits]) << r_off;
  offset += n_bits;

  return code;
}

#endif

#ifdef HAS_ZLIB

int ReadCompressed(gzFile file, char *mem, long start, long num)
{
  // Read data compressed if necessary
  char_type *stackp;
  int finchar;
  code_int code, oldcode, incode;
  char *d,*de;

  // Go to beginning of file
  gzrewind(file);

  // Check whether file is compressed
  if ((gzgetc(file)!='\037') || (gzgetc(file)!=('\235'&0xFF))) {
    // Read data uncompressed
    gzseek(file, start, SEEK_SET);
    gzread(file, mem, num);
    return(1==1);
  }

  // Read data compressed
  d=mem;
  de=mem+num;
  free_ent = 0;			/* first unused entry */
  block_compress = BLOCK_MASK;
  clear_flg = 0;

  maxbits = gzgetc(file);	/* set -b from file */
  block_compress = maxbits & BLOCK_MASK;
  maxbits &= BIT_MASK;
  maxmaxcode = 1L << maxbits;

  if (maxbits > BITS) {
    fprintf(stderr,"compressed with %d bits, can only handle %d bits\n",
            maxbits, BITS);
    return(1==0);
  }

  /*------- initialize the first 256 entries in the table.------*/

  maxcode = MAXCODE(n_bits = INIT_BITS);
  for (code=255;code>=0;code--) {
    codetab[code] = 0;
    tab_suffixof(code) = (char_type)code;
  }
  free_ent = ((block_compress) ? FIRST : 256 );

  finchar = oldcode = getcode(file, 1==1);
  if (oldcode == -1) {/* EOF already? */
    return(1==0);			/* Get out of here */
  }
  *d=finchar;
  start--;
  if (start<0) d++;
  stackp = de_stack;

  while ( ((code = getcode(file, 1==0)) > -1)&&(d<de) ) {
    if ((code == CLEAR)&&block_compress) {
      for (code=255;code>=0;code--)
        codetab[code] = 0;
      clear_flg = 1;
      free_ent = FIRST - 1;
      if ((code=getcode(file, 1==0))==-1) { /* O, untimely death! */
        return(1==0);
      }
    }
    incode = code;

    /*--------- Special case for KwKwK string.------- */

    if (code>=free_ent) {
      *stackp++ = finchar;
      code = oldcode;
    }

    /*------ Generate output characters in reverse order------ */

    while (code>=256) {
      *stackp++ = tab_suffixof(code);
      code = codetab[code];
    }
    *stackp++ = finchar = tab_suffixof(code);

    /*------- And put them out in forward order -------*/
    do {
      *d=*--stackp;
      start--;
      if (start<0) d++;
    }  while (( stackp > de_stack )&&(d<de));

    /*-------- Generate the new entry.------- */

    if ( (code=free_ent) < maxmaxcode ) {
      codetab[code] = (unsigned short)oldcode;
      tab_suffixof(code) = finchar;
      free_ent = code+1;
    }

    /*---------- Remember previous code.------- */

    oldcode = incode;
  }
  return(1==1);
}

#else

int ReadCompressed(FILE *file, char *mem, long start, long num)
{
  // Read data compressed if necessary
  char_type *stackp;
  int finchar;
  code_int code, oldcode, incode;
  char *d,*de;

  // Go to beginning of file
  rewind(file);

  // Check whether file is compressed
  if ((fgetc(file)!='\037') || (fgetc(file)!=('\235'&0xFF))) {
    // Read data uncompressed
    fseek(file, start, SEEK_SET);
    fread(mem, num, 1, file);
    return(1==1);
  }

  // Read data compressed
  d=mem;
  de=mem+num;
  free_ent = 0;			/* first unused entry */
  block_compress = BLOCK_MASK;
  clear_flg = 0;

  maxbits = fgetc(file);	/* set -b from file */
  block_compress = maxbits & BLOCK_MASK;
  maxbits &= BIT_MASK;
  maxmaxcode = 1L << maxbits;

  if (maxbits > BITS) {
    fprintf(stderr,"compressed with %d bits, can only handle %d bits\n",
            maxbits, BITS);
    return(1==0);
  }

  /*------- initialize the first 256 entries in the table.------*/

  maxcode = MAXCODE(n_bits = INIT_BITS);
  for (code=255;code>=0;code--) {
    codetab[code] = 0;
    tab_suffixof(code) = (char_type)code;
  }
  free_ent = ((block_compress) ? FIRST : 256 );

  finchar = oldcode = getcode(file, 1==1);
  if (oldcode == -1) {/* EOF already? */
    return(1==0);			/* Get out of here */
  }
  *d=finchar;
  start--;
  if (start<0) d++;
  stackp = de_stack;

  while ( ((code = getcode(file, 1==0)) > -1)&&(d<de) ) {
    if ((code == CLEAR)&&block_compress) {
      for (code=255;code>=0;code--)
        codetab[code] = 0;
      clear_flg = 1;
      free_ent = FIRST - 1;
      if ((code=getcode(file, 1==0))==-1) { /* O, untimely death! */
        fprintf(stderr,"Can't read compressed file\n");
        return(1==0);
      }
    }
    incode = code;

    /*--------- Special case for KwKwK string.------- */

    if (code>=free_ent) {
      *stackp++ = finchar;
      code = oldcode;
    }

    /*------ Generate output characters in reverse order------ */

    while (code>=256) {
      *stackp++ = tab_suffixof(code);
      code = codetab[code];
    }
    *stackp++ = finchar = tab_suffixof(code);

    /*------- And put them out in forward order -------*/
    do {
      *d=*--stackp;
      start--;
      if (start<0) d++;
    }  while (( stackp > de_stack )&&(d<de));

    /*-------- Generate the new entry.------- */

    if ( (code=free_ent) < maxmaxcode ) {
      codetab[code] = (unsigned short)oldcode;
      tab_suffixof(code) = finchar;
      free_ent = code+1;
    }

    /*---------- Remember previous code.------- */

    oldcode = incode;
  }
  return(1==1);
}

#endif
