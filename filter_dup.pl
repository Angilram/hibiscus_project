while($line = <>){

	if ($line =~ />/){

		$name = $line;
		$HASH{$name} = "";
	}else{

		$HASH{$name} .= $line;
	}
}
foreach $key (keys %HASH){

	print "$key$HASH{$key}";

}
