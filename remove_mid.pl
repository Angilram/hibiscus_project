use Data::Dumper;
#take in fasta
open(file1,$ARGV[0])||die "perl *.pl Fasta.fa Blast.rawblast\n";
while($line = <file1>){

	chomp $line;
	if($line =~ />/){
		$name = ($line =~ m/>(.*)/)[0];
		$HASH{$name} = "";
	}else{
		$HASH{$name} .= $line;
	}
}
#take in blastfile
open(file2,$ARGV[1])||die "perl *.pl Fasta.fa Blast.rawblast\n";
while($line = <file2>){

	chomp $line;
	@array = split "\t", $line;
	$HASH2{$array[2]} = $array[2];

}
#print Dumper(\%HASH2);
foreach $key (keys %HASH2){
	
	print ">$key\n$HASH{$key}\n";

}
