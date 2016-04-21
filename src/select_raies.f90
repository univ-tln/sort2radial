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


subroutine select_raies(mat_inter,tab_dspm_seuil_signal,params)

! The aim of this function is to build the params structure that contains
! the main parameters of the simulation

! Initializations
! ---------------
use init_params_module

implicit none

! Parameters structure
type(parameters)::params
!! Matrice interspectrale sur toutes les raies
! 3eme dimension=toutes les frequences
complex(kind=8),dimension(params%NAT,params%NAT,params%NPTS)::mat_inter
! Intermediate variables to load data
real(kind=8),dimension(params%NPTS)::tab_dspm_seuil_signal,tab_dB_seuil_signal
! Intermediate variables
real(kind=8),dimension(params%NPTS)::dspm,dsp,dspmlog,Lbruit_rapide,dspm_dyn,chameau,dspmconnex,vect_new
complex(kind=8),dimension(params%NPTS)::filtre_BF,dspmlogm,fdspmlog,lp_fdspmlog
real(kind=8),allocatable,dimension(:)::zb,multiseuilref,zbruit_sur_support
integer,allocatable,dimension(:)::sort_index
integer::na,largeurbruit,ierr,ind_noise,i,ind,largref,nzlig,dzlig,Lzero,deltaexca,deltaexcabord
real(kind=8)::pi,snr,noise_level,dB_noise,noiseref
real(kind=8)::dspm_seuil_signal,dspm_seuil_alarm,dB_seuil_signal
real(kind=8)::y1z,y2z,Vexca,y1g,y2g,y1d,y2d,dBmarge
real(kind=8)::dB_dN_bragg,dB_dN_snr,dB_dP_bragg,dB_dP_snr
integer,allocatable,dimension(:)::Izero,Izero_modulo,Ideltag,Ideltad,bfcz
real(kind=8),allocatable,dimension(:)::zero2oneg,zero2oned
integer::x1g,x2g,Ldeltag,x2gg,x2d,Ldeltad,x1dd,x1d
complex::j
integer::lensav,lenwrk,Ldelta,nfc6,nfcz,ind_bruit1,ind_bruit2,lindn,lindp
real(kind=8),allocatable,dimension(:)::work,vect_int
real(kind=8),allocatable,dimension(:)::wsave
real(kind=8)::dB_seuil_dynaN,dB_seuil_dynaP,rap_pic_signal,dB_seuil_picN,dB_seuil_picP,p_seuil_picN,p_seuil_picP,rap_inter_pic
real(kind=8)::qmax,rapanis,NPpicanis,Npicmin,dyna,Nqloc,Ppicmin,Pqloc,dspl_signal_loc,dspm_seuil_dyna
integer::yabossen,yabossep,oklevel,iccloc,igoup,igodown,connex
integer,dimension(1)::icfloc
integer,allocatable,dimension(:)::Nidomaine,Pidomaine

        
! Initialization
pi=acos(0.0)*2.
j=cmplx(0,1)
 
! Calcul du spectre moyen dspm somme sur les antennes
! = somme sur les diagonales
! Boucle sur les antennes
dspm=0.
do na=1,params%NAT
     dsp=abs(mat_inter(na,na,:))
     dspm=dspm+dsp
enddo
! Passage en echelle dB
dspmlog=10.*log10(dspm)

! Seuillage des raies a traiter selon le rapport signal bruit
! Seuillage defini par l'utilisateur params%dB_snr
snr=10.**(params%dB_snr/10.)

! Initialisations diverses    
!Pidomaine=0
!Nidomaine=0

! Test si methode 0 ou 1
if((params%methselect==0) .OR. (params%methselect==1))then

     ! Seuil 1 : dB_noise
     ! !!!!!!!!!!!!!!!!!!
     ! Nombre de raies considerees pour determiner le niveau de bruit
     largeurbruit=size(params%liste_raies_bruit)
     allocate(zb(largeurbruit),stat=ierr)
     allocate(sort_index(largeurbruit),stat=ierr)
     allocate(zbruit_sur_support(largeurbruit),stat=ierr)
     
     ! Niveau du bruit apres classement croissant du dspm hors des bandes de Bragg
     ! et hors de la bande centrale
     zb=dspm(params%liste_raies_bruit)
     sort_index=(/(i, i=1,largeurbruit, 1)/)
     call ssort(zb,sort_index,largeurbruit)

     ! Premier seuil possible: la valeur moyenne de zb
     noise_level=0.
     do ind_noise=1,largeurbruit
         noise_level=noise_level+zb(ind_noise)
     enddo
     noise_level=noise_level/largeurbruit
     dB_noise=10.*log10(noise_level)
     
     ! Seuil 2 : dspm_seuil_alarm
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!   
     ! Autre calcul de seuil par definition probabiliste du seuil par
     ! minimum des prolongements de 3 lois "erf" calees en 3 points
     ! differents de l histogramme zb (lois gaussiennes pour la puissance et
     ! non pour l'amplitude). Akaike?
     allocate(multiseuilref(size(params%proba_ref)),stat=ierr)

     ! Boucle sur les valeurs de depart
     do ind=1,size(params%proba_ref)
         ! Indice de depart
         largref=floor(largeurbruit*params%proba_ref (ind))
    
         ! Valeur correspondante
         noiseref=zb(largref)
    
         ! Prolongement probabiliste gaussien à 0.95!
!         kmargin=scipy.special.erfinv(1-params%pfalarm)/scipy.special.erfinv(params%proba_ref[ind])
         
         ! Stockage
         multiseuilref(ind)=params%kmargin(ind)*noiseref
     enddo

     ! On garde le seuil le plus bas
     dspm_seuil_alarm=minval(multiseuilref)

     ! Seuil 3 : dB_seuil_signal
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!   
     ! La methode statistique du dspm_seuil_alarm est plus subtile.
     ! On prend pour seuil, le minimum des deux seuils 1 et 2.mod_select_raies.py
     ! snr*noise_level = rapport signal sur bruit defini par l'utilisateur x
     ! niveau de bruit moyen sur le signal.
     dspm_seuil_signal=min(snr*noise_level,dspm_seuil_alarm)

     ! Echelle log
     dB_seuil_signal=10.*log10(dspm_seuil_signal)

     ! Seuil 4 : tab_dspm_seuil_signal
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Seuil variable selon allure spectre
     ! Fabrication en 3 temps d'un socle courbe seuil lisse sous la courbe log
     dspmlogm=dspmlog

     ! Socle lineaire double largeur sous la raie zero large
     nzlig=nint(params%kzero_wide*params%NPTS)
     dzlig=nint(nzlig/2.)
     allocate(Izero(2*dzlig+1),stat=ierr)
     allocate(Izero_modulo(2*dzlig+1),stat=ierr)
     Izero=(/(i, i=-dzlig,dzlig, 1)/)
     Lzero=size(Izero)
     Izero_modulo=modulo(Izero,params%NPTS)+1
     y1z=dspmlogm(Izero_modulo(1))
     y2z=dspmlogm(Izero_modulo(2*dzlig+1))
     ! Interpolation lineaire pour remplacer la zone proche du zero
     ! Doppler
     dspmlogm(Izero_modulo)=y1z+((y2z-y1z)*(Izero+dzlig)/(Lzero-1))

     ! Pour passer largement sous les regions de Bragg elargies
     ! borne liee à vmax elargie
     Vexca=max(2., (1.4*params%vmax))
     deltaexca=nint((2.*Vexca)/(params%Lambda*params%DeltaF))
     deltaexcabord=deltaexca
     ! Test pour ne pas depasser le 0 doppler
     if (deltaexcabord>(params%DIbragg-1))then
         deltaexcabord=params%DIbragg-1
     endif
     
     ! Socle a gauche lineaire sous la bosse de Bragg 
     x1g=((params%NPTS+1)-params%DIbragg-deltaexca)
     x2g=((params%NPTS+1)-params%DIbragg+deltaexcabord)
     ! Test au cas où
     if (x2g>params%NPTS)then
         x2g=params%NPTS
     endif

     ! Socle_gauche=x1g:x2g;
     allocate(Ideltag(deltaexca+deltaexcabord+1),stat=ierr)
     allocate(zero2oneg(deltaexca+deltaexcabord+1),stat=ierr)
     Ideltag=(/(i, i=-deltaexca,deltaexcabord, 1)/)
     Ldeltag=size(Ideltag)
     zero2oneg=(Ideltag+deltaexca)/(Ldeltag-1.)
     ! Petite moyenne autour des raies limites
     y1g=0.
     do ind=x1g-5,x1g
         y1g=y1g+dspmlogm(ind)
     enddo
     y1g=y1g/6.
     
     x2gg=x2g+5
     if(x2gg>params%NPTS)then
         x2gg=params%NPTS
     endif
     
     y2g=0.
     do ind=x2g,x2gg
         y2g=y2g+dspmlogm(ind)
     enddo
     y2g=y2g/(x2gg-x2g+1)

     ! Mise à jour du dspmlogm
     dspmlogm((params%NPTS+1)-params%DIbragg+Ideltag)=y1g+((y2g-y1g)*zero2oneg)

     ! Socle a droite lineaire sous la bosse de Bragg
     x1d=1+params%DIbragg-deltaexcabord
     if (x1d<1)then
         x1d=1
     endif
     
     x2d=1+params%DIbragg+deltaexca
     ! Socle_droit=x1d:x2d;
     allocate(Ideltad(deltaexca+deltaexcabord+1),stat=ierr)
     allocate(zero2oned(deltaexca+deltaexcabord+1),stat=ierr)
     Ideltad=(/(i, i=-deltaexcabord,deltaexca, 1)/)
     Ldeltad=size(Ideltad)
     zero2oned=(Ideltad+deltaexcabord)/(Ldeltad-1.)
     
     ! Petite moyenne autour des raies limites
     x1dd=x1d-5
     if (x1dd<1)then
         x1dd=1
     endif
     
     y1d=0
     do ind=x1dd,x1d
         y1d=y1d+dspmlogm(ind)
     enddo
     y1d=y1d/(x1d-x1dd+1)
     
     y2d=0.
     do ind=x2d,x2d+5
         y2d=y2d+dspmlogm(ind)
     enddo
     y2d=y2d/6.
     
     ! Mise à jour du dspmlogm
     dspmlogm(1+params%DIbragg+Ideltad)=y1d+((y2d-y1d)*zero2oned)

     ! Lissage de la courbe seuil par filtrage passe bas FFT 
     ! (sinon segments de droites + bosses eventuelles bateaux...):
     ! Spectre
     lensav=2*params%NPTS+int(log(real(params%NPTS,kind=8)))+4
     lenwrk=2*params%NPTS
     allocate(wsave(lensav),stat=ierr)
     allocate(work(lenwrk),stat=ierr)
     call zfft1i(params%NPTS,wsave,lensav,ierr)
     
     fdspmlog=dspmlogm
     
     call zfft1f(params%NPTS,1,fdspmlog,params%NPTS,wsave,lensav,work,lenwrk,ierr)         

     fdspmlog=fdspmlog*params%NPTS

     ! Choix de la frequence passe bas de coupure -6dB.
     ! La courbe ne doit pas suivre les bosses de Bragg
     ! Attention il y en a autant de part et d autre
     ! prenons la largeur de l excavation Ldelta
     Ldelta=nint((Ldeltag+Ldeltad)/2.)
     nfc6=floor(params%NPTS/real(Ldelta,kind=8))

     ! Frequence d'attenuation totale
     nfcz=2*nfc6

     ! Fabrication du filtre a phase nulle
     allocate(bfcz(nfcz),stat=ierr)
     bfcz=(/(i, i=1,nfcz, 1)/)

     filtre_BF(bfcz)=0.5*(1.+cos(pi*(bfcz-1)/nfcz))
     filtre_BF(params%NPTS+1-bfcz)=0.5*(1.+cos(pi*(bfcz)/nfcz))

     ! Application du filtre
     lp_fdspmlog=filtre_BF*fdspmlog

     ! Calcul de la courbe de seuil lisse qui forme un socle
     call zfft1b(params%NPTS,1,lp_fdspmlog,params%NPTS,wsave,lensav,work,lenwrk,ierr)
     
     tab_dB_seuil_signal=real(lp_fdspmlog,kind=8)/params%NPTS

     ! Evaluation subtile du bruit rapide au dessus du lisse
     ! On pourrait faire un calcul avec Pfa et Pref
     Lbruit_rapide=abs(dspmlogm-tab_dB_seuil_signal)

     ! Tri sur la zone pur bruit     
     zbruit_sur_support=Lbruit_rapide(params%liste_raies_bruit)
     sort_index=(/(i, i=1,largeurbruit, 1)/)
     call ssort(zbruit_sur_support,sort_index,largeurbruit)
     ind_bruit1=floor(largeurbruit/3.)
     ind_bruit2=floor(2.*largeurbruit/3.)
     
     ! dBmarge=prolongement asymptotyque de la distribution et accentuation
     ! Evolution possible : tenir compte du nombre d echantillons
     ! independants
     dBmarge=((2.+params%accent)*zbruit_sur_support(ind_bruit2)-zbruit_sur_support(ind_bruit1))
     ! dBmarge compris entre params%dB_snr et 2*params%dB_snr
     dBmarge=max(dBmarge,params%dB_snr)
     dBmarge=min(dBmarge,2.*params%dB_snr)

     ! On rajoute la marge de bruit
     tab_dB_seuil_signal=tab_dB_seuil_signal+dBmarge
     ! et dBmarge est devenu automatique

     ! Fabrication du seuil en echelle puissance lineaire
     tab_dspm_seuil_signal=10.**(tab_dB_seuil_signal/10.)

     ! Synthese des 3 tests : prend-on a gauche ou a droite la version de
     ! seuil variable?
     ! Test a gauche
     ! Valeur max du signal sur la zone de Bragg
     dB_dN_bragg=maxval(dspmlog(params%NPTS-params%DIbragg+Ideltag))
     ! Valeur du signal le plus faible en principe 
     dB_dN_snr=dB_dN_bragg-dB_noise
     ! Test sur la dynamique maximale totale entre la plus forte raie 
     ! et le plus faible signal accepte pour le domaine.
     ! Si le rapport signal sur bruit est bon on garde les seuils statiques
     if(dB_dN_snr>params%dB_dynaraie)then
         ! Application du seuil 3
         allocate(vect_int(params%NPTS-params%raie_neg_min+1),stat=ierr)
         vect_int=1.
         tab_dspm_seuil_signal((/(i, i=params%raie_neg_min,params%NPTS, 1)/))=dspm_seuil_signal*vect_int
         ! Echelle dB    
         tab_dB_seuil_signal((/(i, i=params%raie_neg_min,params%NPTS, 1)/))=dB_seuil_signal*vect_int
         deallocate(vect_int)
     endif

     ! Test a droite
     ! Valeur max du signal sur la zone de Bragg
     dB_dP_bragg=maxval(dspmlog(params%DIbragg+Ideltad))
     ! Valeur du signal le plus faible en principe 
     dB_dP_snr=dB_dP_bragg - dB_noise
     ! Test sur la dynamique maximale totale entre la plus forte raie 
     ! et le plus faible signal accepte pour le domaine.
     ! Si le rapport signal sur bruit est bon on garde les seuils statiques
     if(dB_dP_snr>params%dB_dynaraie)then
         ! Application du seuil 3       
         allocate(vect_int(params%raie_pos_max),stat=ierr)
         vect_int=1.
         tab_dspm_seuil_signal((/(i, i=1,params%raie_pos_max, 1)/))=dspm_seuil_signal*vect_int
         ! Echelle dB    
         tab_dB_seuil_signal((/(i, i=1,params%raie_pos_max, 1)/))=dB_seuil_signal*vect_int
     endif

     ! Definition par defaut des seuils de dynamique
     dB_seuil_dynaN=tab_dB_seuil_signal(params%NPTS-params%DIbragg)
     dB_seuil_dynaP=tab_dB_seuil_signal(params%DIbragg)

     ! Seuillage des pics, centre de bosses spectrales de Bragg, 
     ! selon le rapport pic sur seuil_signal.
     ! On est plus exigeant que pour une raie de signal de base
     ! puisqu a partir de ce point on selectionne un voisinage.
     ! L'objectif est d'eviter les remontees vers les contributions du 2eme
     ! ordre (contraste entre 1.5dB et 2.5dB)
     rap_pic_signal=10.**(params%dB_new_pic/10.)

     ! Nouveaux seuils dynamiques
     dB_seuil_picN=dB_seuil_dynaN+params%dB_new_pic
     dB_seuil_picP=dB_seuil_dynaP+params%dB_new_pic

     ! Puissance seuil pour acceptation d'un nouveau pic 
     ! premiere initialisation
     p_seuil_picN=10.**(dB_seuil_picN/10.)
     p_seuil_picP=10.**(dB_seuil_picP/10.)

     ! Dynamique maximale entre les pics locaux du meme cote (sens des vagues)
     ! Idee au depart pour discriminer l ordre 2.
     ! Critere + lache maintenant
     rap_inter_pic=10.**(params%dB_inter_pic/10.)

     ! Boucle de recherche des bosses
     ! ------------------------------------------
     ! Initialisation du spectre Dynamique de recherche qui sera ensuite excave    
     dspm_dyn=dspm

     ! On attenue ce qui n est pas dans la zone selectionnee support_vague
     ! hors_jeu=find(params%support_vague==0);
     ! On attribue une valeur faible
     do ind=1,params%NPTS
         if(params%support_vague(ind)==0.)then
             dspm_dyn(ind)=1.
         endif
     enddo

     ! Recherche du pic maximum maximorum dans dspm_dyn initial complet
     qmax=maxval(dspm_dyn)

     ! Anisotropie max droite gauche (entre raie doppler negative et
     ! positive)
     rapanis=10.**(params%dB_anis/10.)
 
     ! Seuil pour admettre la recherche de raie sur le cote faible
     NPpicanis=qmax/rapanis

     ! Recherche des pics de bosses et de leurs voisinages
     ! D'abord recherche a gauche
     ! Compteur de bosses du cote negatif
     yabossen=0         
     ! Variable test de sortie de niveau des pics locaux
     oklevel=1

     ! Nouveau seuil minimum pour declarer une nouvelle bosse
     Npicmin=max(NPpicanis,p_seuil_picN)

     ! Dynamique maximale du premier ordre
     dyna=10.**(params%dB_dynaraie/10.)
     
     vect_new=0

     ! Boucle sur les bosses
     do while((yabossen<params%nbbosse) .AND. (oklevel==1))
         ! Copie du spectre excavable
         chameau=dspm_dyn
    
         ! Recherche du plus grand max de ce cote et de sa position
         Nqloc=maxval(chameau((/(i, i=params%raie_neg_min,params%NPTS, 1)/)))
         icfloc=maxloc(chameau((/(i, i=params%raie_neg_min,params%NPTS, 1)/)))

         ! Position dans le grand tableau
         iccloc=icfloc(1)+params%raie_pos_max
    
         ! Puissance seuil pour acceptation d'un nouveau pic    
         p_seuil_picN=tab_dspm_seuil_signal(iccloc)*rap_pic_signal
    
         ! Nouveau seuil minimum pour declarer une nouvelle bosse
         Npicmin=max(Npicmin,p_seuil_picN)  

         ! Verification du niveau du pic et sortie sinon        
         if(Nqloc<Npicmin)then
             ! Il n'y a plus rien a prendre de ce cote
             oklevel=0
         else
             ! Seuillage dynamique du lobe du premier ordre    
             dspl_signal_loc=Nqloc/dyna
             dB_seuil_dynaN=10.*log10(dspl_signal_loc)
        
             ! Si premiere bosse
             if(yabossen==0)then
                 ! Nouvelle initialisation du seuil lateral de pic
                 Npicmin=max(Nqloc/rap_inter_pic,p_seuil_picN)
             endif
             
             ! Selection de la bosse intervalle locale
             dspmconnex=1.

             ! On demarre de la position du max locale
             ! On va vers les indices croissants
             igoup=iccloc
        
             ! Indice sur la tolerance de connexite
             connex=params%toler_con
        
             ! On cherche la limite sup de zone de bosse
             do while((connex>=1) .AND. (igoup<=params%NPTS))
                 dspm_seuil_dyna=max(tab_dspm_seuil_signal(igoup),dspl_signal_loc)
                 ! Test si sous le seuil et limite de connexite
                 if (dspm_dyn(igoup)<=dspm_seuil_dyna)then 
                     connex=connex-1 
                     dspmconnex(igoup)=1
                 else
                     dspmconnex(igoup)=dspm_dyn(igoup)
                 endif
                 
                 ! Le point igoup est traite
                 dspm_dyn(igoup)=1
                 igoup=igoup+1 
             enddo
             
             !  On va vers les indices decroissants
             igodown=iccloc-1
        
             ! Indice sur la tolerance de connexite
             connex=params%toler_con
        
             ! On cherche la limite basse de la zone de bosse
             do while((connex>=1) .AND. (igodown>=params%raie_neg_min))
                 dspm_seuil_dyna=max(tab_dspm_seuil_signal(igodown),dspl_signal_loc)
                 ! Test si sous le seuil et limite de connexite
                 if (dspm_dyn(igodown)<=dspm_seuil_dyna)then
                     connex=connex-1 
                     dspmconnex(igodown)=1
                 else
                     dspmconnex(igodown)=dspm_dyn(igodown)
                 endif
                 
                 ! Le point igodown est traite
                 dspm_dyn(igodown)=1
                 igodown=igodown-1
             enddo
             
             ! Mise a jour du domaine
             do ind=1,params%NPTS
                 if(dspmconnex(ind)>1)then
                     vect_new(ind)=1.
                 endif
             enddo
        
             ! Increment du nombre de bosses
             yabossen=yabossen+1
         endif
     enddo
     
     lindn=0

     ! Tri sur le domaine negatif
     do ind=1,params%NPTS
         if(vect_new(ind)>0)then
             lindn=lindn+1
         endif
     enddo
     
     allocate(Nidomaine(lindn),stat=ierr)
     lindn=0
     do ind=1,params%NPTS
         if(vect_new(ind)>0)then
             lindn=lindn+1
             Nidomaine(lindn)=ind
         endif
     enddo

     ! Recherche sur le cote positif
     ! --------------------------------------
     ! Compteur de bosses
     yabossep=0 

     ! Nouveau seuil minimum pour declarer une nouvelle bosse
     Ppicmin=max(NPpicanis,p_seuil_picP)

     ! Variable de sortie d'apres le test de niveau des pics locaux
     oklevel=1

     vect_new=0

     ! Boucle sur les bosses
     do while ((yabossep<params%nbbosse) .AND. (oklevel==1))
         ! Copie du spectre excavable
         chameau=dspm_dyn

         ! Recherche du plus grand max de ce cote et de sa position
         Pqloc=maxval(chameau((/(i, i=1,params%raie_pos_max, 1)/)))
         icfloc=maxloc(chameau((/(i, i=1,params%raie_pos_max, 1)/)))
         iccloc=icfloc(1)
    
         ! Puissance seuil pour acceptation d'un nouveau pic    
         p_seuil_picP=tab_dspm_seuil_signal(iccloc)*rap_pic_signal
    
         ! Nouveau seuil minimum pour declarer une nouvelle bosse
         Ppicmin=max(Ppicmin,p_seuil_picP)
         
         ! Verification du niveau du pic et sortie sinon
         if(Pqloc<Ppicmin)then
             oklevel=0
         else
             ! Seuillage dynamique du lobe du premier ordre    
             dspl_signal_loc=Pqloc/dyna
             dB_seuil_dynaP=10.*log10(dspl_signal_loc)        

             ! Test si premiere bosse
             if (yabossep==0)then
                 ! Nouveau seuil minimum pour declarer une nouvelle
                 ! bosse
                 Ppicmin=max(Pqloc/rap_inter_pic,p_seuil_picP)
             endif
            
             ! Selection de la bosse locale
             dspmconnex=1.
             
             ! On demarre de la position du max locale
             ! On va vers les indices croissants
             igoup=iccloc
        
             ! Indice sur la tolerance de connexite
             connex=params%toler_con
        
             ! On cherche la limite sup de la zone de bosse
             do while ((connex>=1) .AND. (igoup<=params%raie_pos_max))
                 dspm_seuil_dyna=max(tab_dspm_seuil_signal(igoup),dspl_signal_loc)
                 ! Test si sous le seuil et limite de connexite
                 if(dspm_dyn(igoup)<=dspm_seuil_dyna)then
                     connex=connex-1
                     dspmconnex(igoup)=1.
                 else
                     dspmconnex(igoup)=dspm_dyn(igoup)
                 endif
                 
                 ! Le point igoup est traite
                 dspm_dyn(igoup)=1.
                 igoup=igoup+1 
             enddo
             
             ! On demarre de la position du max locale
             ! On va vers les indices decroissants
             igodown=iccloc-1
        
             ! Indice sur la tolerance de connexite
             connex=params%toler_con
        
             ! On cherche la limite basse de la zone de bosse
             do while((connex>=1).AND.(igodown>=0))
                 dspm_seuil_dyna=max(tab_dspm_seuil_signal(igodown),dspl_signal_loc)
                 ! Test si sous le seuil et limite de connexite
                 if (dspm_dyn(igodown)<=dspm_seuil_dyna)then
                     connex=connex-1 
                     dspmconnex(igodown)=1.
                 else
                     dspmconnex(igodown)=dspm_dyn(igodown)
                 endif
                 
                 ! Le point igodown est traite
                 dspm_dyn(igodown)=1.
                 igodown=igodown-1     
             enddo
             
             ! Mise a jour du domaine
             do ind=1,params%NPTS
                 if(dspmconnex(ind)>1)then
                     vect_new(ind)=1.
                 endif
             enddo        
             
             ! Increment du nombre de bosses
             yabossep=yabossep+1
             
         endif
     enddo    
     
     lindp=0
      ! Tri sur le domaine negatif
     do ind=1,params%NPTS
         if(vect_new(ind)>0)then
             lindp=lindp+1
         endif
     enddo
     
     allocate(Pidomaine(lindp),stat=ierr)
     lindp=0
     do ind=1,params%NPTS
         if(vect_new(ind)>0)then
             lindp=lindp+1
             Pidomaine(lindp)=ind
         endif
     enddo
 
     allocate(params%idomaine(lindn+lindp),stat=ierr)
     
     do ind=1,lindp
         params%idomaine(ind)=Pidomaine(ind)
     enddo
     
     do ind=1,lindn
         params%idomaine(lindp+ind)=Nidomaine(ind)
     enddo
      
endif
 
return
end subroutine select_raies
