#include <stdio.h>
#include <stdlib.h>

	 struct node
	{ int key; struct node *next; };
	static struct node *head, *z, *t; 
/*        struct node *head, *z, *t; */

#ifdef _UNDERSCORE
void stackinit_ ()
#else
void	stackinit() 
#endif
	   {
	     head = (struct node *) malloc(sizeof *head);
	     z = (struct node *) malloc(sizeof *z); 

/*             head = (struct node *) malloc(sizeof(head));
             z = (struct node *) malloc(sizeof(z)); */


             head->next = z; head->key=0;
	     z->next = z;
	     z->key = 0;
	   }

#ifndef _UNDERSCORE
void	push(p)
#else
void push_ (p)
#endif
           int *p;
	   {
	     int v;
	     v = *p;
	     t = (struct node *) malloc(sizeof *t);	
	     t->key = v; t->next = head->next;	
	     head->next =t;	
	   }

#ifndef _UNDERSCORE
void	pop(x)
#else
void pop_ (x)
#endif
           int *x;
	   {
	     t = head->next; head->next = t->next;
	     *x = t->key;
	     free(t);
	   }

#ifndef _UNDERSCORE
void	stackempty(i)
#else
void stackempty_ (i)
#endif
          int *i;
	  { 
	    *i = 0;
            if(head->next == z) *i = 1;
          }

#ifndef _UNDERSCORE
void        stackflush()
#else
void stackflush_ ()
#endif
           {
             free(head);
             free(z);
           }
