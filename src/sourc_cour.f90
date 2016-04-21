! Copyright Actimar 2016

! This software called SORT2RADIAL is a computer program whose purpose is to 
! derive radial current velocities from input SORT or USORT Wera radar 
! data using the MUSIC algorithm. It is a fortran version of 
! the Matlab code developped by "Laboratoire de Sondages Electromagnétiques
! de l'Environnement Terrestre, UMR 6017" now integrated into the MIO
! (Mediterranean Institute of Oceanography).

! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software.  You can  use, 
! modify and/ or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info". 

! As a counterpart to the access to the source code and  rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty  and the software's author,  the holder of the
! economic rights,  and the successive licensors  have only  limited
! liability.

! In this respect, the user's attention is drawn to the risks associated
! with loading,  using,  modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean  that it is complicated to manipulate,  and  that  also
! therefore means  that it is reserved for developers  and  experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or 
! data to be ensured and,  more generally, to use and operate it in the 
! same conditions as regards security.

! The fact that you are presently reading this means that you have had
! knowledge of the CeCILL license and that you accept its terms.


subroutine sourc_cour(data_chevfft,var_ant,norm_vp_noise,n_porte,dopplertp_m,params)

! The aim of this function is to build the params structure that contains
! the main parameters of the simulation

! Initializations
! ---------------
use init_params_module

implicit none

! Parameters structure
type(parameters)::params
! Data after Doppler processing
complex(kind=8),dimension(params%NPTS,params%n_ant,params%n_OBS,params%n_range)::data_chevfft
! Antenna manifold
complex(kind=8),dimension(params%NAT,params%nazim)::var_ant
! Intermediate variables to load data
real(kind=8),dimension(params%prof,params%n_dir,params%multiv2f,2)::dopplertp_m
real(kind=8),dimension(params%NAT)::norm_vp_noise,vp_noise6,diag_gammas
integer::n_porte
! Matrice interspectrale sur toutes les raies
! 3eme dimension=toutes les frequences
complex(kind=8),allocatable,dimension(:,:,:)::mat_inter
complex(kind=8),dimension(params%NAT,params%NAT)::gamma,mgamma,gammaS,mat_int
! Vecteurs observations (a une frequence donnee),
! par antenne, et pour chaque porte et chaque
! echantillon d observation
complex(kind=8),allocatable,dimension(:,:,:)::fobs
complex(kind=8),allocatable,dimension(:,:)::obs,subspace_select,projvsurs
real(kind=8),dimension(params%n_dir,params%multiv2f,2)::doppler
integer,allocatable,dimension(:)::inter_liss
real(kind=8),dimension(params%NPTS)::tab_dspm_seuil_signal,tab_sort_dspm
real(kind=8),dimension(params%NPTS,params%NAT)::ligne_vpip
integer,dimension(params%NAT+1,params%NPTS)::jazX
real(kind=8),dimension(params%NAT+1,params%NPTS)::jnivX
real(kind=8),dimension(params%NAT,params%NAT)::xr,xi,mgammar,mgammai
complex(kind=8),dimension(params%NAT,params%NAT)::vecproX
real(kind=8),dimension(params%NAT)::wr,wi
complex,dimension(params%NAT)::valpro
real(kind=8),dimension(params%NAT)::valproX,dB_valproX
integer,dimension(params%NPTS)::vect_int,pos1
real(kind=8),dimension(params%nazim)::s2col,denom,qf,ql,qlmax,qlmax0,qlmax_up,sort_ql
real(kind=8),allocatable,dimension(:)::lnorm_vp_noise,nrjSOURCE
integer,dimension(3)::indajust
logical::solution_libre
! Intermediate variables
integer::pbase,n_echant,n_dim,n_space,ierr,nf,i_r,j_r,i,ind_noise,&
     nlsur2,n_liss,n_seek,indice_niveau_ql,n_max,nb_source,compteur_horschamp,i_s
integer::iporte,ibase,jbase,i_ant,nbf,nb_nomax,total_crown,nid,nbajust,ivptest,&
     nb_source_0,n_seek_0,ind_ref,ind,compteur_source,i_freq,ns,jazim,j11,mv1
real(kind=8)::pi,chouilla_diag,noise_level,snr_collage,leval_sigma,eval_sigma,&
     sigma2moyen,denommin,maxteta,tetas,dtetas
real(kind=8),allocatable,dimension(:)::ko_seek_max,qf_int,nrjSOURCE_test,&
     ampSOURCE_test,ko_mean_score,ko_9_score
integer,allocatable,dimension(:)::interseek,imax,vect_test,i_musicmax,&
     position_angulaire_test,intersource,solution_angulaire
real(kind=8),allocatable,dimension(:,:)::solution_angulaire_seek,&
     nrjSOURCE_seek,lnrjSOURCE_seek,proj2,ACOVSOURCE_test,ko_rap7,alco_ref
complex(kind=8),allocatable,dimension(:,:)::A,pinvA,COVSOURCE_test,&
     intA,matPROD_test,vect_intl,vect_intc
real(kind=8)::ideal_beam_width,beam_width_site,w_pied_source,niveau_ql,&
     seuil_ql,ko_max,ko_9
complex::j
real(kind=8)::ref_mean,ref_max
integer,dimension(1)::nb_source_pos

! Initialization
pi = acos(0.0)*2.
j=cmplx(0,1)

allocate(mat_inter(params%NAT,params%NAT,params%NPTS),stat=ierr)

! Indices importants
! Indice de la premiere porte incluse de la couronne
pbase=n_porte-floor(params%reg_radial/2.)

! Nombre d'observations echantillons qui sont sommes en interspectral
n_echant=params%reg_radial*params%n_OBS
n_dim   =params%NAT

! Specification de la derniere vp non nulle pour les futures cadrages
n_space=min(n_dim,params%q_OBS)

! Vecteurs observations (a une frequence donnee),
! par antenne, et pour chaque porte et chaque
! echantillon d observation
allocate(fobs(params%NAT,n_echant,params%NPTS),stat=ierr)
allocate(obs(params%NAT,n_echant),stat=ierr)

! Preparation du tableau Doppler resultat pour la couronne (porte)
! avec la resolution angulaire de la grille cherchee doppler multivalue
! Variable qui sert en fin de fonction a ecrire dans la variable globale
! dopplertp = doppler "toutes portes"
doppler=999999999.

chouilla_diag=0.001    

! Premiere boucle sur toutes les frequences
! Calcul des matrices interspectrales pour toutes les frequences doppler
do nf=1,params%NPTS
     ! Construction des echantillons observes
     ! pbase est retabli comme premiere porte incluse dans le groupe
     do i_r=1,params%reg_radial
         iporte=(pbase-1)+i_r
         ibase=(i_r-1)*params%n_OBS        
         do j_r=1,params%n_OBS
             jbase=ibase+j_r
             obs(:,jbase)=data_chevfft(nf,params%listant,j_r,iporte)
         enddo
     enddo
     
     ! On en aura besoin plus tard pour la formation de voie vers les
     ! sources!
     fobs(:,:,nf)=obs

     ! Calcul effectif de gamma, la matrice interspectrale pour la frequence
     ! Pour eviter les problemes de convergence on va mettre un petit
     ! qqchose sur la diagonale.
     gamma=0.
     do i_ant=1,params%NAT
         gamma(i_ant,i_ant)=chouilla_diag
     enddo

     ! Somme des echantillons via produit matriciel
     gamma=gamma+matmul(obs, transpose(conjg(obs))) 

     !normalisation pour garder des amplitudes comparables
     gamma=gamma/(params%NAT*n_echant)
     mat_inter(:,:,nf)=gamma
enddo

! Selection de raies pour le traitement DF (idomaine) avec retour de
! tab_dspm_seuil_signal pour l'evaluation du bruit
call select_raies(mat_inter,tab_dspm_seuil_signal,params)

! Definition de noise level pour la couronne
! Tri sur les niveaux de seuil
vect_int=(/(i, i=1,params%NPTS, 1)/)
call ssort(tab_dspm_seuil_signal,vect_int,params%NPTS)
tab_sort_dspm=tab_dspm_seuil_signal

! Position empirique pour determiner le bruit (!)
ind_noise=ceiling(size(tab_dspm_seuil_signal)/32.)
noise_level=tab_sort_dspm(ind_noise)

! Valeur pour le collage des vp de bruit sur les vp calculees.
! Avec une bonne trajectoire des vp on peut se rapprocher!
! Le /2 est empirique, le bruit sur la trajectoire des vp de bruit est lie
! au nombre d echantillons de meme que le bruit sur les spectres
snr_collage=10.**((params%dB_snr/2.)/10.)

! Nombre de raies selectionnees
nbf=size(params%idomaine)

! Matrice locale de toutes valeurs propres mise a un au depart
ligne_vpip=1.

! Direction multi source (incremente de un pour stabilite boucle)
jazX=999999999

! Niveau multi source (incremente de un pour stabilite boucle)
jnivX=999999999.

! Accumulateur du nombre de facteur music sans max valides
nb_nomax=0

! Variable intermediaire pour compter le nombre de sources par couronne
total_crown=0

! Initialisation pour le lissage du facteur music
nlsur2=floor(params%n_liss0/2.)
allocate(inter_liss(2*nlsur2+1),stat=ierr)
inter_liss=(/(i, i=-nlsur2,nlsur2, 1)/)
n_liss=size(inter_liss)

! Seconde boucle sur les raies detectees.
! Diagonalisation seulement sur le domaine selectionne et recherche
! directionnelle des sources.
do nid=1,nbf !nbf
     ! Raie correspondante
     nf=params%idomaine(nid)

     ! Matrice interspectrale de la raie
     mgamma=mat_inter(:,:,nf)

     ! Diagonalisation de mgamma
     mgammar=real(mgamma,kind=8)
     mgammai=aimag(mgamma)
     
     call cg(params%NAT,params%NAT,mgammar,mgammai,wr,wi,1,xr,xi,ierr)
     
     valpro=wr+j*wi
     vecproX=xr+j*xi
     
     ! Mise en format colonne des valeurs propres
     valproX=abs(valpro)
     
     ! On remplit la ligne de la matrice memoire
     ligne_vpip(nf,:)=valproX

     ! Colonne echelle dB des valeurs propres
     dB_valproX=10.*log10(valproX)

     ! Determination de la dimension de l espace source valide 
     ! (ns=nombre de sources) pour le calcul ensuite du facteur MUSIC.
     ! On normalise vp_noise sur les 3 petites valeurs propres
      ! qui sont le plus souvent (toujours) dans l espace bruit
      ! exemple: mn_space=8, on regarde sur les 3 plus petites valeurs
      ! puis comparaison a la courbe comme dans 4 mais avec le decalage.
 
      ! Seuil "qui marche" permet de differencier cas avec beaucoup
      ! d'echantillons et/ou antennes.
     allocate(lnorm_vp_noise(n_space),stat=ierr)
     if(n_space>=6)then
         lnorm_vp_noise=10.*log10(norm_vp_noise(1:n_space))
         ! Ajustement trajectoire bruit modelisee sur vp reelles (en dB)
         ! Dans ce cas on travaille sur les 3 dernieres
         indajust=(/(i, i=n_space-2,n_space, 1)/)
         nbajust=size(indajust)
         leval_sigma=sum(dB_valproX(indajust))-sum(lnorm_vp_noise(indajust))
         leval_sigma=leval_sigma/nbajust
         eval_sigma=10.**(leval_sigma/10.)
     else
         !sinon sur la derniere : a ameliorer : adapter intervale 
         eval_sigma=valproX(n_space)/norm_vp_noise(n_space)
     endif
     
     ! Calcul trajectoire seuil finale pour detection sources
     vp_noise6=snr_collage*eval_sigma*norm_vp_noise

     ! Initialisation
     nb_source_0=0
     ivptest=1

     ! Recherche nombre sources : on doit etre < mnspace
     ! On s'arrete quand la valeur des valeurs propres reelles devient
     ! inferieure a celle du bruit modelisee
     do while((valproX(ivptest)>vp_noise6(ivptest)) .AND. (ivptest<n_space))              
         nb_source_0=nb_source_0+1
         ivptest=ivptest+1                               
     enddo
     
     ! Test inutile en general
     if(nb_source_0>=params%NAT)then
         nb_source_0=params%NAT-1
     endif

     ! Pour garder au moins une solution si snr bon
     if(valproX(1)>(2*noise_level))then
         nb_source_0=max(nb_source_0,1)
     endif

     ! Majoration du nombre de sources considerees ! Surement a revoir!
     ! n_seek_0 est un maximum. On va chercher des sources de 1 a n_seek_0.
     n_seek_0=min(nb_source_0,params%nb_source_max)
     n_seek_0=max(n_seek_0,1)

     ! Initialisation des tableaux pour test comparatif
     allocate(ko_mean_score(n_seek_0),stat=ierr)
     allocate(ko_9_score(n_seek_0),stat=ierr)
     allocate(ko_seek_max(n_seek_0),stat=ierr)
     allocate(solution_angulaire_seek(n_seek_0,n_seek_0),stat=ierr)
     allocate(nrjSOURCE_seek(n_seek_0,n_seek_0),stat=ierr)
     allocate(lnrjSOURCE_seek(n_seek_0,n_seek_0),stat=ierr)

     ko_seek_max=1.
     solution_angulaire_seek=999999999.
     nrjSOURCE_seek=999999999.
     lnrjSOURCE_seek=999999999.

     ! Boucle sur le nombre croissant n_seek de sources recherchees
     do n_seek=1, n_seek_0
         ! Domaine
!         allocate(interseek(n_seek),stat=ierr)
!         interseek=(/(i, i=1,n_seek, 1)/)
!         interseek=range(2)
         allocate(interseek(2),stat=ierr)
         interseek=(/(i, i=1,2, 1)/)
         
         ! Calcul de sigma2, niveau moyen du bruit sur les valeurs
         ! propres faibles qui ne correspondent pas a des sources.
         sigma2moyen=sum(valproX((/(i, i=3,params%NAT, 1)/)))/params%NAT

         ! Calcul du facteur music qf
         ! ----------------------------------
         ! Selection des vecteurs propres qui nous interessent
         allocate(subspace_select(params%NAT,2),stat=ierr)
         subspace_select=vecproX(:,interseek)
         
         ! On projette la variete (reponse normalisee) sur le sous
         ! espace. x projete sur y se calcule par y'x.
         allocate(projvsurs(2,params%nazim),stat=ierr)
         projvsurs=matmul(transpose(conjg(subspace_select)), var_ant)

         ! On prend la somme des carres des projections.
         ! Mise a jour du denominateur sans post normalisation qui a
         ! deja ete faite dans rep_ant.
         allocate(proj2(2,params%nazim),stat=ierr)
         proj2=conjg(projvsurs)*projvsurs
         s2col=sum(proj2,1)
         
         !         interseek=range(n_seek)
         deallocate(interseek)
         allocate(interseek(n_seek),stat=ierr)
         interseek=(/(i, i=1,n_seek, 1)/)
         
         ! Calcul de denominateur
         denom=1.
         denom=denom-s2col

         ! Pre-test de validite du resultat du calcul
         denommin=minval(denom)
         if (denommin<=0)then
             ! Modification de denom. Impact?
             denom=denom-denommin+0.001
         endif

         ! Facteur music qf
         qf=1./denom
         
!         open(ind_file_file,file=trim('qf.bin'),access='direct',recl=8*720)
!         write(ind_file_file,rec=1)qf
!         close(ind_file_file) 

         ! Validation et evaluation de la solution a n_seek sources 
         ! Recherche "subtile" des max de fmusic selon la convexite de
         ! fmusicX.
         ! Determination des "pseudomax" locaux de fmusicX.    
         ! Quand la bosse est large ou si deux bosses se recouvrent
         ! il n y a qu un max sur ql
         ! On va donc chercher les max sur (music - music lisse)

         ! Calcul de la version lisse du facteur music  
         ! Lissage circulaire : yes please!
         call smooth(qf,inter_liss,n_liss,params%nazim,ql)
!         open(ind_file_file,file=trim('ql.bin'),access='direct',recl=8*720)
!         write(ind_file_file,rec=1)ql
!         close(ind_file_file) 
         ! Recherche des max par convexite
         ! convexite locale = brut -lisse
         qlmax0=qf-ql

         ! Recherche des max locaux de qlmax0 en deux passes
         call pass1(qlmax0,params%nazim,qlmax_up)
         call pass2(qlmax_up,params%nazim,qlmax)

         ! Suppression des max symetriques des reseaux alignes
         if (params%ya_sym)then
             qlmax(int(params%hors_jeu_sym))=0.
         endif
         
         ! Calcul du seuil de validite des max music
         ! Determination du nombre d'elements concernes par les raies de
         ! Bragg en fonction de la largeur du faisceau empirique
         ideal_beam_width=2./params%n_ant
         beam_width_site=params%beam_width*ideal_beam_width
         w_pied_source=params%nazim/(2.*pi)*beam_width_site

         ! Concerne le facteur music lisse
         ! Tri croissant
         sort_ql=ql
         allocate(vect_test(params%nazim),stat=ierr)
         vect_test=(/(i, i=1,params%nazim, 1)/)
         call ssort(sort_ql,vect_test,params%nazim)
         deallocate(vect_test)

         indice_niveau_ql=params%nazim-(n_seek*nint(w_pied_source))
         ! Forcement sur la partie superieure sinon le seuil descend trop
         ! bas.
         indice_niveau_ql=max((params%nazim/2),indice_niveau_ql)
         niveau_ql=sort_ql(indice_niveau_ql)
         ! Extrapolation par la pente de la zone hors source
         seuil_ql=niveau_ql*(params%nazim/(indice_niveau_ql+1.))   

         ind_ref=0
        
         ! Identification et comptage des max valide (>seuil_ql)
         do ind=1,params%nazim
             if((qlmax(ind)>0.).AND.(qf(ind)>seuil_ql))then
                 ind_ref=ind_ref+1
             endif
         enddo
    
         allocate(imax(ind_ref),stat=ierr)
         ind_ref=0
         do ind=1,params%nazim
             if((qlmax(ind)>0.).AND.(qf(ind)>seuil_ql))then
                 ind_ref=ind_ref+1
                 imax(ind_ref)=ind
             endif
         enddo

         if(ind_ref==0)then
             n_max=0
             nb_nomax=nb_nomax+1                  
         else
             n_max=ind_ref
         endif

         ! Si le nombre de max trouves est inferieur au nombre de
         ! sources cherchees
         ! A revoir : n_max<n_seek cette solution devrait egalement etre evaluee?
         if (n_max>=n_seek)then
             ! On travaille seulement sur les n_seek plus grands max du
             ! facteur music.

             ! Classement decroissant des "pseudomax" de music.
             ! Le classement reste aussi utile pour le cas a 1 source
             allocate(qf_int(n_max),stat=ierr)
             allocate(i_musicmax(n_max),stat=ierr)
             allocate(position_angulaire_test(n_seek),stat=ierr)
             
             qf_int=qf(imax)
             i_musicmax=(/(i, i=1,n_max, 1)/)
             call ssort_dec(qf_int,i_musicmax,n_max)
             position_angulaire_test=imax(i_musicmax(interseek))
             deallocate(i_musicmax)
             deallocate(qf_int)
            
             ! Variete d'antenne pour la position angulaire trouvee pour les
             ! raies
             allocate(A(params%NAT,n_seek),stat=ierr)
             allocate(intA(n_seek,n_seek),stat=ierr)
             allocate(pinvA(n_seek,params%NAT),stat=ierr)
             
             A=var_ant(:,position_angulaire_test)
             
             ! Matrice pseudo inverse de A
             call cinv(matmul(transpose(conjg(A)), A),n_seek,n_seek,intA)
             pinvA=matmul(intA, transpose(conjg(A)))

             ! Facteur de correction de l estimation de l energie du
             ! bruit par les vp faibles.
             do ind=1,params%NAT
                 mat_int(ind,ind)=sigma2moyen
             enddo
             gammaS=mgamma-mat_int       
            
             ! On verifie que gammaS n'a pas d elements diagonaux
             ! negatifs et sinon on les rend positif petits
             do ind=1,params%NAT
                 diag_gammas(ind)=real(gammaS(ind,ind),kind=8)
                 if(diag_gammas(ind)<0.)then
                     gammaS(ind,ind)=0.01
                 endif
             enddo
             
             allocate(COVSOURCE_test(n_seek,n_seek),stat=ierr)
             allocate(ACOVSOURCE_test(n_seek,n_seek),stat=ierr)
             allocate(matPROD_test(n_seek,n_seek),stat=ierr)
             allocate(ko_rap7(n_seek,n_seek),stat=ierr)
             allocate(nrjSOURCE_test(n_seek),stat=ierr)
             allocate(ampSOURCE_test(n_seek),stat=ierr)
             allocate(vect_intl(1,n_seek),stat=ierr)
             allocate(vect_intc(n_seek,1),stat=ierr)
             
             ! Matrice covariance des sources
             COVSOURCE_test=matmul(matmul(pinvA, gammaS), transpose(conjg(pinvA)))

             ! Estimation d independace
             ACOVSOURCE_test=abs(COVSOURCE_test)
             
             do ind=1,n_seek
                 nrjSOURCE_test(ind)=ACOVSOURCE_test(ind,ind)
             enddo
             
             ! Analyse de covariance des sources seulement si n_seek>1
             if(n_seek>1)then
                 ! Construction de la matrice seuil acceptation
                 ampSOURCE_test=sqrt(nrjSOURCE_test)
                 do ind=1,n_seek
                     vect_intl(1, ind)=ampSOURCE_test(ind)
                     vect_intc(ind, 1)=ampSOURCE_test(ind)
                 enddo
                 matPROD_test=matmul(vect_intc, conjg(vect_intl))

                 ! Pour evaluation statistique de ko_test
                 ! Normalisation des covariances par rapport a l'amplitude
                 ! des covariances
                 ko_rap7=ACOVSOURCE_test/matPROD_test
                
                 ! Determination de la valeur max de covariance
                 ko_max=maxval(ko_rap7)
                
                 ! Stockage pour test ulterieur
                 ko_seek_max(n_seek)=ko_max
             endif
             
             ! Mise en memoire de la solution courante dans tableaux
             solution_angulaire_seek(n_seek,interseek)=position_angulaire_test
             nrjSOURCE_seek(n_seek,interseek)=nrjSOURCE_test
             lnrjSOURCE_seek(n_seek,interseek)=10.*log10(nrjSOURCE_test)

             deallocate(position_angulaire_test)
             deallocate(A)
             deallocate(pinvA)
             deallocate(intA)
             deallocate(COVSOURCE_test)
             deallocate(ACOVSOURCE_test)
             deallocate(nrjSOURCE_test)
             deallocate(matPROD_test)
             deallocate(ko_rap7)
             deallocate(ampSOURCE_test)
             deallocate(vect_intl)
             deallocate(vect_intc)             
         endif
         deallocate(imax)
         deallocate(interseek)
         deallocate(subspace_select)
         deallocate(proj2)
         deallocate(projvsurs)

         ! Test si le nombre de sources cherchees est > 1
         if (n_seek_0>1)then
             ! Test de la meilleure valeur du tableau ko_seek_max
             ! 1/ Construction tableau complet alcov_ref interpole selon
             ! q_OBS.
             ! Initialisation : nombre de sources max cherchees x (ko_mean,ko_max)
             allocate(alco_ref(n_seek_0,2),stat=ierr)
             alco_ref=1.
             
             ! Boucle sur le nombre de sources recherchees
             do i_s=2, n_seek_0
                 call read_alcov2(i_s,ref_mean,ref_max,params)
                 alco_ref(i_s,1)=ref_mean
                 alco_ref(i_s,2)=ref_max
             enddo
             
             ! 2/ Evaluation du score de ko_seek_max par rapport a ko_mean
             ! et ko_max.
             ko_mean_score=ko_seek_max/alco_ref(:,1)
             ko_9_score=ko_seek_max/alco_ref(:,2)
    
             ! 3/ Meilleur score sur ko_mean
             nb_source_pos=minloc(ko_mean_score)
             nb_source=nb_source_pos(1)
             
             ! 4/ Choix par rapport a source unique
             ! Si toutes les solutions sont mauvaises (score vs 9ieme decile >1)
             ! => on ne garde qu'une source
             ko_9=ko_9_score(nb_source)
             if(ko_9>1)then
                 nb_source=1
             endif
             deallocate(alco_ref)
         else
             nb_source=1
         endif
         
         ! Une fois nb_source choisi, ca roule:
         allocate(intersource(nb_source),stat=ierr)
         allocate(solution_angulaire(nb_source),stat=ierr)
         allocate(nrjSOURCE(nb_source),stat=ierr)
         intersource=(/(i, i=1,nb_source, 1)/)
         solution_angulaire=solution_angulaire_seek(nb_source,intersource)
         nrjSOURCE=nrjSOURCE_seek(nb_source,intersource)

         ! Calcul de l'energie si le nombre de sources est >0 et remplissage de
         ! tableaux
         if (nb_source>0)then   
             ! Pour le moment pas de choix dans la methode de calcul de l'energie.
             ! On fait confiance au calcul sur la matrice de covariance
!             PsourceX=nrjSOURCE
            
             ! Enregistrement des resultats de la raie dans les tableaux
             ! couronnes
             !print nid, nf, intersource
             jazX(intersource,nf)=solution_angulaire
             jnivX(intersource,nf)=nrjSOURCE
             total_crown=total_crown+nb_source
         endif
         deallocate(intersource)
         deallocate(solution_angulaire)
         deallocate(nrjSOURCE)
     enddo
     deallocate(lnorm_vp_noise)
     deallocate(ko_seek_max)
     deallocate(ko_mean_score)
     deallocate(ko_9_score)
     deallocate(solution_angulaire_seek)
     deallocate(nrjSOURCE_seek)
     deallocate(lnrjSOURCE_seek)
enddo

! Statistiques des sources
if (nbf>=1)then
     params%total_raies=params%total_raies+nbf
     params%total_sources=params%total_sources+total_crown
endif

! Mise en forme des solutions dans le tableau local "doppler"
maxteta=modulo((params%tetad(params%n_dir)-params%tetad(1)),360.)
compteur_source=0
compteur_horschamp=0
pos1=999999999

! Traitement quand il y a un domaine selectionne non vide (nombre de raies
! Doppler detectees >0)
if (nbf>0)then
     ! Boucle sur les raies Doppler detectees
     do i_freq=1,nbf
         ! Indice de la raie Doppler detectee
         nf=params%idomaine(i_freq)
         ns=1
         ! Remplissage de doppler
         do while ((jazX(ns,nf)<999999999) .AND. (ns<=params%NAT))
             ! Recuperation de l'indice azimutal
             jazim=jazX(ns,nf)
        
             ! Valeur de l'azimut en degre
             tetas=params%azimd(jazim)
        
             ! Azimut dans le repere du radar
             dtetas=modulo(tetas-params%tetad(1),360.)
        
             if ((dtetas>=0) .AND. (dtetas<=maxteta))then
                 ! Il est plus correct de prendre l arrondi
                 ! On ajoute 1 car les indices commencent a 1
                 ! Indice de l'azimut dans le repere du radar
                 j11=nint(dtetas/params%dteta)+1
                 if (ns==1)then
                     pos1(nf)=params%tetad(j11)
                 endif
                 
                 ! Recherche de l'indice libre sur le tableau doppler pour
                 ! remplir a la suite les valeurs
                 mv1=1
                 solution_libre=(mv1<=params%multiv2f).AND.(doppler(j11,mv1,1)==999999999.)
                 do while ((.NOT.solution_libre).AND.(mv1<params%multiv2f))
                     mv1=mv1+1
                     solution_libre=(mv1<=params%multiv2f).AND.(doppler(j11,mv1,1)==999999999.)
                 enddo
                 if (solution_libre)then
                     doppler(j11,mv1,1)=nf
                     ! Voila la ligne a changer
                     ! Ne pas attribuer utiliser directement la  relation croissante
                     ! fmusic et valeurs propres
                     ! Donc utiliser le calcul de niveau de la source
                     doppler(j11,mv1,2)=jnivX(ns,nf)                   
                     compteur_source=compteur_source+1
                 else
!                     disp(' ')
!                     disp('DEPASSEMENT de CAPACITE SOURCES D''UNE CELLULE')
!                     disp('Changer params%multiv2f pour ce site!!!!!')
                 endif
             else !if dtetas
                 compteur_horschamp=compteur_horschamp+1
             endif
             ns=ns+1   
         enddo            
     enddo
endif

! Mise a jour du tableau pour toutes les couronnes
do ind=1,params%prof
     if(params%porte(ind)==n_porte)then
         ind_ref=ind
         exit
     endif
enddo
dopplertp_m(ind_ref,:,:,:)=doppler(:,:,:)

deallocate(params%idomaine)
deallocate(mat_inter)
deallocate(fobs)
deallocate(obs)
deallocate(inter_liss)

return
end subroutine sourc_cour
