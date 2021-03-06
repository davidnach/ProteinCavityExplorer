/* name: utility.c              */
/* author: J. Michael Word      */
/* date written: 2/26/96        */
/* purpose: utility functions   */

/*****************************************************************/
/* NOTICE: This is free software and the source code is freely   */
/* available. You are free to redistribute or modify under the   */
/* conditions that (1) this notice is not removed or modified    */
/* in any way and (2) any modified versions of the program are   */
/* also available for free.                                      */
/*               ** Absolutely no Warranty **                    */
/* Copyright (C) 1999 J. Michael Word                            */
/*****************************************************************/

#include <stdio.h>
#include <ctype.h>

void note(char *message) {
	fprintf(stderr, "%s\n", message);
}

void warn(char *message) {
	fprintf(stderr, "WARNING: %s\n", message);
}

void errmsg(char *message) {
	fprintf(stderr, "ERROR: %s\n", message);
}

void halt(char *message) {
	fprintf(stderr, "ERROR: %s\n", message);
	exit(1);
}

int parseInteger(char *str, int start, int len) {
	register int value = 0;
	register char ch;
	int neg = 0, inside = 0;

	if (!str || start < 0) { return 0; }
	str += start;

	while((len-- > 0) && *str) {
		ch = *str++;
		if ((ch >='0') && (ch <= '9')) {
			value = (10*value) + (ch - '0');
			inside = 1;
		}
		else if (ch == '+' && !inside) {
			inside = 1;
		}
		else if (ch == '-' && !inside) {
			neg = 1;
			inside = 1;
		}
		else if (isspace(ch) && !inside) { /* nothing */ }
		else break; /* end of integer */
	}
	return (neg?-value:value);
}

float parseReal(char *str, int start, int len) {
	double value = 0.0, scale = 1.0;
	register char ch;
	int inside = 0, infract = 0;

	if (!str || start < 0) { return 0; }
	str += start;

	while((len-- > 0) && *str) {
		ch = *str++;
		if ((ch >='0') && (ch <= '9')) {
			value = (10.0*value) + (ch - '0');
			if (infract) { scale *= 0.1; }
			inside = 1;
		}
		else if (ch == '+' && !inside) {
			inside = 1;
		}
		else if (ch == '-' && !inside) {
			scale = -1.0;
			inside = 1;
		}
		else if (ch == '.' && !infract) {
			infract = 1;
		}
		else if (isspace(ch) && !inside) { /* nothing */ }
		else break; /* end of real */
	}
	return value*scale;
}

void copyChars(char *to, char *from, int n) {
	int i;

	for(i=0; i<n; i++) { to[i] = from[i]; }
}

int compArgStr(char *str, char *arg, int min) {
	int i, max;
	char s, a;

	if (!str || !arg) return 0;

	max = strlen(arg);

	for(i=0; i<max; i++) {
		s = toupper(str[i]);
		a = toupper(arg[i]);

		if (i >= min && (s == '\0' || s == '.'
			|| s == '+' || s == '-' || isdigit(s))) {
			break; /* good ending point */
		}
		else if (s != a) {
			i = 0; /* failed to match */
			break;
		}
	}

	return i;
}
