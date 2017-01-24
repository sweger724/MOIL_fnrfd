	subroutine mdiv_add_les(totmon,mdivlist,n_cpy)
	include 'COMMON/LENGTH.BLOCK'
	integer totmon,mdivlist(0:maxdiv+1,maxmono)
	integer n_cpy(*)

c local
	integer i,j,k,i_save
	integer mdivlist_tmp(0:maxdiv+1,maxmono)

	mdivlist_tmp(1,1) = 0
        do 290 i=1,totmon
                do  291 j=1,mdivlist(0,i)
                        i_save = mdivlist_tmp(j,i)
                        do 292 k=mdivlist(j,i)+1,mdivlist(j+1,i)
                        i_save = i_save + n_cpy(k)
292             continue
                mdivlist_tmp(j+1,i) = i_save
291             continue
                mdivlist_tmp(1,i+1) = mdivlist_tmp(mdivlist(0,i)+1,i)
290     continue

        do 293 i=1,totmon
                do 293 j=1,mdivlist(0,i)+1
                        mdivlist(j,i) = mdivlist_tmp(j,i)
293     continue

	return
	end
