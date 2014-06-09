#include <config.h>
#include<stdio.h>

int main()
{ 
  double db=0.0;
  int it=0;
  unsigned int uit=0;
  string text;
  float ft=0.0;

  config tmp("../config/testconfig.cfg");


    cout << "euro " << tmp.getParameter("euro",db) << endl; 

    cout << "eindrittel " << tmp.getParameter("eindrittel",ft) << endl; 

    cout << "pi " << tmp.getParameter("pi",ft) << endl;

    cout << "tausend " << tmp.getParameter("tausend",it) << endl;

    if(tmp.existParameter("text")) {
      cout << "Text " << tmp.getParameter("text",text) << endl; }

    cout << "Integernum " << tmp.getParameter("integernum",it) << endl;

    cout << "Unsigned int " << tmp.getParameter("integernum",uit) << endl;

    cout << "minuseins " << tmp.getParameter("minuseins",it)  << endl;


  it=-3;
  if(tmp.setParameter("minuseins",it)){
    cout << "minuseins was updated with setParameter" << endl;}

  cout << "minuseins " << tmp.getParameter("minuseins",it)  << endl;

  uit=42;

  if(tmp.changeParameter("integernum",uit)){
    cout << "integernum was changed !" << endl; }

  if(tmp.delParameter("text")) {
    cout << "text wurde geloescht" << endl;}



 if(tmp.delParameter("floatzahl")) {
    cout << "floatzahl wurde geloescht" << endl; }



   if(tmp.newParameter("testconfig1.cfg","testnew",ft))
    {
      cout << "The parameter testnew was created" << endl;
      cout << "testnew " << tmp.getParameter("testnew",ft) << endl;
    }
 
   ft=0;
   string txt="neuer text :) 1.12.99";
   if(tmp.newParameter("testconfig.cfg","neutext",txt)) {
     cout << "neutext lautet " << tmp.getParameter("neutext",txt) << endl; };
     return 0;
}






