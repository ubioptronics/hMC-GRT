#ifndef __TRMC_OPTS_H__
#define __TRMC_OPTS_H__
/***************************************************
							Global #defines; 
***************************************************/
#define MCSILENT (fopen("/dev/null", "w") == NULL) ? stdout : fopen("/dev/null", "w")
#define MCDEBUG MCSILENT
#define die(s) { perror(s); exit(EXIT_FAILURE); }
#define MALLOC(s,t) if(((s) = malloc(t)) == NULL) { die("error: malloc() "); }
#define DELIMITERS " ,:;/|"
//***************************************************/
#define MUA 0
#define MUS 1
#define ANIS_G 2
#define REFIND 3
#define MUA_F 0
#define TAU 1
#define FLQY 2

#endif
