#define _GNU_SOURCE
#define __USE_GNU
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <byteswap.h>
#include <string.h>

#ifndef DATADECODE
#define DATADECODE fa250DataDecode
#endif

char progName[255];
void usage();

int
main(int argc, char *argv[])
{
  int iarg = 1;
  int byteswap = 0;
  FILE *f;
  uint32_t data = 0;
  char *user_input = NULL;
  size_t input_size = 0;
  size_t line_len = 0;
  int bytes_read, bytes_skip = 0;


  /* Evaluate the command line arguments */
  strncpy((char *) &progName, argv[0], 255);

  if(argc > 1)
    {
      /* help */
      if ((strcmp(argv[iarg],"-h") == 0) || (strcmp(argv[iarg],"--help") == 0))
	{
	  usage();
	  return 0;
	}

      /* byteswapping specified */
      if ((strcmp(argv[iarg],"-b") == 0) || (strcmp(argv[iarg],"--byteswap") == 0))
	{
	  byteswap = 1;
	  ++iarg;
	}

      /* Filename included */
      if (iarg == argc-1)
	{
	  f = fopen(argv[iarg], "r");
	  if(f)
	    {
	      printf("\n");
	      while(fscanf(f, "%x", &data) > 0)
		{
		  if(byteswap)
		    data = bswap_32(data);

		  DATADECODE(data);
		}
	    }
	  else
	    perror(argv[1]);

	  return 0;
	}
    }

  /* User data entry */
  printf("\n");
  while(1)
    {
      line_len = getline(&user_input, &input_size, stdin);
      if(line_len > 1)
	{
	  while(sscanf(user_input + bytes_skip, "%x%n", &data,
		       &bytes_read) == 1)
	    {
	      if(byteswap)
		data = bswap_32(data);

	      DATADECODE(data);
	      bytes_skip += bytes_read;
	    }
	  bytes_skip = 0;
	}
      else
	break;
    }

  return 0;
}

void
usage()
{
  printf("\nUsage: ");
  printf("%s [OPTION] [FILENAME]\n", progName);
  printf("If [FILENAME] is not specified, user input is accepted.\n");
  printf("\n");
  printf("OPTIONs:\n");
  printf("  -h, --help                    This help page.\n");
  printf("  -b, --byteswap                Byteswap (32bit) before decoding.\n");
  printf("\n\n");


}
