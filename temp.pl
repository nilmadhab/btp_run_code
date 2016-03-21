for($w=15; $w <=50 ; $w=$w+5)
{
	print $w;
	print "------------------";
	system("perl generate.pl $w 140 0");
	
}
