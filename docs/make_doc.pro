pro make_doc, outdir=outdir
   
   ;; build HTML documentation  an output directory can be
   ;; specified with OUTDIR
    
    idir = getenv('HSRED_DIR')

    spawn, 'find '+idir+'/idl/spec2d -name "*.pro"  -print', prolist1
    spawn, 'find '+idir+'/idl/spec1d -name "*.pro"  -print', prolist2
    spawn, 'find '+idir+'/idl/preproc -name "*.pro"  -print', prolist3
    spawn, 'find '+idir+'/idl/utils -name "*.pro"  -print', prolist4
    

    prolist = [prolist1,prolist2,prolist3,prolist4]
    
    if n_elements(outdir) eq 0L then outdir = idir+'/docs/'
    
    make_html_help, prolist, outdir+'hsred_doc.html', $
      /strict, /link_files, title='IDL Help for HSRED', /verbose

return
end
