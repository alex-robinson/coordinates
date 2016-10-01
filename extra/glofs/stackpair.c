#include <stdio.h>
#include <stdlib.h>

	 struct node
	{ int key; struct node *next; };
	 struct node *heada, *za, *ta;
	 struct node *headb, *zb, *tb;

#ifndef _UNDERSCORE
void	stackpairinit() 
#else
void stackpairinit_ () 
#endif
	   {
	     heada = (struct node *) malloc(sizeof *heada);
	     za = (struct node *) malloc(sizeof *za);
	     heada->next = za; heada->key=0;
	     za->next = za;
	     za->key = 0;

	     headb = (struct node *) malloc(sizeof *headb);
	     zb = (struct node *) malloc(sizeof *zb);
	     headb->next = zb; headb->key=0;
	     zb->next = zb;
	     zb->key = 0;
	   }

#ifndef _UNDERSCORE
void	stackpairflush() 
#else
void stackpairflush_ () 
#endif
	   {
             free(heada);
             free(headb);
             free(za);
             free(zb);
	   }

#ifndef _UNDERSCORE
void	pushpair(pa, pb)
#else
void pushpair_ (pa, pb)
#endif
             int *pa;
             int *pb;
	   {
	     int va;
	     int vb;
	     va = *pa;
	     ta = (struct node *) malloc(sizeof *ta);	
	     ta->key = va; ta->next = heada->next;	
	     heada->next =ta;	

	     vb = *pb;
	     tb = (struct node *) malloc(sizeof *tb);	
	     tb->key = vb; tb->next = headb->next;	
	     headb->next =tb;	
	   }

#ifndef _UNDERSCORE
void	poppair(xa,xb)
#else
void poppair_ (xa,xb)
#endif
             int *xa;
             int *xb;
	   {
	     ta = heada->next; heada->next = ta->next;
	     *xa = ta->key;
	     free(ta);
	     tb = headb->next; headb->next = tb->next;
	     *xb = tb->key;
	     free(tb);
	   }

#ifndef _UNDERSCORE
void	stackpairempty(i)
#else
void stackpairempty_ (i)
#endif
            int *i;
	  { 
	    *i = 0;
            if(heada->next == za) *i = 1;
          }
