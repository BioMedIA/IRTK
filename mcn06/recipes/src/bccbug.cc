/* The following function is recommended by Borland Technical Support to
   "fix" the error "Floating Point Formats Not Linked".  To use this file,
   compile it along with your own files on the compiler command line. You
   do not need to call it, just compile it along with your files.         */

void LinkFloat(void)
{
  float a=0, *b=&a;
  a=*b;
}
