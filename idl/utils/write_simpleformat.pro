PRO write_simpleformat, rootname

    sphectfile = findfile('spHect-' + rootname + '.fits*')
    hs_toiraf, sphectfile[0], 'spIraf-' + rootname+'.fits'

    spzbestfile = findfile('spZbest-' + rootname + '.fits*')
    hs_spZbest_to_ascii, spzbestfile[0], 'spZbest-' + rootname + '.dat'

END
