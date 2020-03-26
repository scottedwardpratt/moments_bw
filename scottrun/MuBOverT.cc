double GetMuBOverT(double roots){
	double a,b,p1,p2,x1,x2,muBoverT;
	x1=1.0/9.1;
	x2=1.0/200;
	p1=log(0.0075);
	p2=log(0.78);
	b=((p1/x1)-(p2/x2))/(x1-x2);
	a=(p1-b*x1*x1)/x1;
	//printf("a=%g, b=%g\n",a,b);
	muBoverT=-0.5*((a/roots)+(b/(roots*roots)));
	/*  // for checking fit
	double poverpbar_STAR[8]={0.769,0.708,0.469,0.058,0.028,0.0078,0.0038,0.0013};
	double roots_STAR[8]={200,130,62.4,17.27,12.32,8.76,7.61,6.27};
	double poverpbar;
	for(int i=0;i<8;i++){
		roots=roots_STAR[i];
		poverpbar=exp(((a/roots)+(b/(roots*roots))));
		printf("roots=%g, muBoverT_STAR=%g =? %g\n",roots,poverpbar_STAR[i],poverpbar);
	}
	*/
	return muBoverT;
}
