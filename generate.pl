$w = $ARGV[0];
$t = $ARGV[1];
$h = $ARGV[2];
#for($w=5; $w <=50 ; $w=$w+5)
{
	system("g++ file_read_write.cpp -o a.out");
	system("./a.out ".$w." ".$t." 5 ".$h);
	#system("g++ 3d_scheduling_dvfs.cpp -o scale");
	#$test = `./a.out 1`;
	#print $test;
	#print "nil";
	system("g++ 3d_scheduling_dvfs.cpp -o dvfs");
	#$test = `./a.out 1`;
	#print $test;
	#print "nil";
	system("./dvfs 1 15000 >". $w."_15000_dvfs_1.txt ");
	#system("./scale 1 25000 >". $w."_25000.txt ");
	
	#system("./scale 1 15000 >". $w."_15000_dvfs.txt ");
	#system("./scale 1 25000 >". $w."_25000.txt ");
}
