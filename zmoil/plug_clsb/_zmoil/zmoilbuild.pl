sub clsbDataCopy {
	my( $dstDir ) = @_;
	
	# COPY
	mkdir( "$dstDir/_zmoil" );
	mkdir( "$dstDir/docs" );
	my @glob = <$zlabDir/../plug_clsb/_zmoil/docs/*.pdf>;
	foreach( @glob ) {
		cp( $_, "$dstDir/docs" );
	}

	my @glob = <$zlabDir/../plug_clsb/_zmoil/*.zui>;
	foreach( @glob ) {
		cp( $_, "$dstDir/_zmoil" );
	}
}

true;
