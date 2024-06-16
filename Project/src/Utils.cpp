#include"Utils.hpp"
#include <cassert>

using namespace Eigen;

namespace LibraryDFN{



//***********************************************PARTE 1***********************************************

//funzione che suddivide la frattura in triangoli (assumo che i vertici del poligono non siano allineati)
std::vector<std::array<unsigned int,3>> triangola_frattura(DFN& disc_frac_net,unsigned int idr) //idr=id relativo, nel senso che è l'indice del vettore numVertici
{
    std::vector<std::array<unsigned int,3>> triangolazione;
    triangolazione.reserve(disc_frac_net.numVertici[idr]-2);
    for(unsigned int i=1;i<disc_frac_net.numVertici[idr]-1;i++)
    {
        std::array<unsigned int,3> triangolo={0,i,i+1};
        triangolazione.push_back(triangolo);
    }
    return triangolazione;
}




//funzione che trova il versore perpendicolare al piano che contiene il poligono
inline Vector3d versore_normale(const std::vector<Vector3d>& poligono)
{
    Vector3d versore(0.,0.,0.);
    for(unsigned int i=1;i<poligono.size()-1;i++)
    {
        versore+=1./(poligono.size()-2)*(poligono[i]-poligono[0]).cross(poligono[i+1]-poligono[0]);
    }
    versore.normalize();
    return versore;
}
//funzione che trova il versore perpendicolare al piano che contiene il poligono, pos=posizione nel vettore DFN.vertici della frattura
inline Vector3d versore_normale_2(const DFN& disc_frac_net,const std::vector<unsigned int>& poligono,const unsigned int& pos)
{
    Vector3d versore(0.,0.,0.);
    for(unsigned int i=1;i<poligono.size()-1;i++)
    {
        versore+=1./(poligono.size()-2)*(disc_frac_net.vertici[pos][poligono[i]]-disc_frac_net.vertici[pos][poligono[0]]).cross(disc_frac_net.vertici[pos][poligono[i+1]]-disc_frac_net.vertici[pos][poligono[0]]);
    }
    versore.normalize();
    return versore;
}




//funzione che trova se un punto è interno o sul bordo del poligono
std::array<bool,2> interno_bordo_poligono(const DFN& disc_frac_net, const std::vector<unsigned int>& poligono,const Vector3d& x,const unsigned int& pos_frattura, const double& tol)
{
    Vector3d n=versore_normale_2(disc_frac_net,poligono,pos_frattura);
    double indicatore=((x-disc_frac_net.vertici[pos_frattura][poligono[poligono.size()-1]]).cross(disc_frac_net.vertici[pos_frattura][poligono[0]]-disc_frac_net.vertici[pos_frattura][poligono[poligono.size()-1]])).dot(n);
    bool interno=(indicatore >-tol);
    bool bordo=false;
    unsigned int i=0;
    while(interno && i<poligono.size())
    {
        indicatore=((x-disc_frac_net.vertici[pos_frattura][poligono[i]]).cross(disc_frac_net.vertici[pos_frattura][poligono[i+1]]-disc_frac_net.vertici[pos_frattura][poligono[i]])).dot(n);
        if(!(abs(indicatore)>tol))
        {
            bordo=true;
        }
    }
    std::array<bool,2> risultato={interno,bordo};
    return risultato;
}




//funzione che scarta le fratture che sicuramente non si intersecano, l'output è un vettore di indici degli indici di DFN.idfratture che non sono stati scartati, quindi su cui bisogna controllare manualmente se si intersechino
std::vector<std::array<unsigned int,2>> scarta_fratture(DFN& disc_frac_net,const double& tol)
{
    std::vector<std::array<unsigned int,2>> indici_validi; //vettore finale di output con gli indici che non hanno passato i test

    //NOTA IMPORTANTE: si potrebbe più facilmente riservare al vettore indici_validi (numFratture su 2) posti, tuttavia se numFratture è troppo alto, potrei andare in overflow e avere un comportamento indefinito
    std::list<std::array<unsigned int,2>> temp; //lista temporanea che contiene le coppie di indici che poi vanno in indici_validi
    std::vector<std::tuple<unsigned int,Vector3d,double>> info_baricentro; //vettore con tutte le info da confrontare per effettuare i test
    info_baricentro.reserve(disc_frac_net.numFratture);

    for(unsigned int i=0;i<disc_frac_net.numFratture;i++)
    {
        //trovo il baricentro della frattura
        Vector3d baricentro(0,0,0);
        for(Vector3d& coord:disc_frac_net.vertici[i])
        {
            baricentro+=coord;
        }
        baricentro=baricentro/disc_frac_net.numVertici[i];
        //trovo il raggio (al quadrato) della sfera che circonscrive il poligono
        double raggio=0.;
        for(Vector3d& coord:disc_frac_net.vertici[i])
        {
            double candidato_raggio=(coord-baricentro).squaredNorm();
            if(candidato_raggio>raggio)
                raggio=candidato_raggio;
        }
        std::tuple<unsigned int,Vector3d,double> tupla(i,baricentro,raggio);
        info_baricentro.push_back(tupla);
    }

    //i test effettivi
    //ciclo su tutte le coppie possibili di fratture
    for(unsigned int i=0;i<disc_frac_net.numFratture-1;i++)
    {
        for(unsigned int j=i+1;j<disc_frac_net.numFratture;j++)
        {
            //test 1: vedo se le sfere che circoscrivono il poligono si intersecano...
            if((std::get<Vector3d>(info_baricentro[i])-std::get<Vector3d>(info_baricentro[j])).squaredNorm()<std::get<double>(info_baricentro[i])+std::get<double>(info_baricentro[j]))
            {
                //...se si intersecano, vedo se riesco a trovare una retta dove applicare il SAT
                double val_min=INFINITY;
                double val_max=-INFINITY;
                for(auto& vertice:disc_frac_net.vertici[i])
                {
                    double proiezione=(vertice-std::get<Vector3d>(info_baricentro[j])).dot(disc_frac_net.versori[j]);
                    val_min=std::min(val_min,proiezione);
                    val_max=std::max(val_max,proiezione);
                }
                // se non vale l'ipotesi del SAT per il primo asse, provo il secondo asse
                if(!(val_min > -tol || val_max < tol))
                {
                    val_min=INFINITY;
                    val_max=-INFINITY;
                    for(auto& vertice:disc_frac_net.vertici[j])
                    {
                        double proiezione=(vertice-std::get<Vector3d>(info_baricentro[i])).dot(disc_frac_net.versori[i]);
                        val_min=std::min(val_min,proiezione);
                        val_max=std::max(val_max,proiezione);
                    }
                    //se non vale l'ipotesi del SAT neanche per il secondo asse, allora c'è una possibile intersezione tra le 2 fratture
                    if(!(val_min > -tol || val_max < tol))
                    {
                        std::array<unsigned int,2> coppia={i,j};
                        temp.push_back(coppia);
                    }
                }
            }
        }
    }
    //sposto tutte le coppie dalla lista al vettore
    indici_validi.reserve(temp.size());
    for(auto& coppia:temp)
    {
        indici_validi.push_back(coppia);
    }
    return indici_validi;
}




//funzione che genera la matrice di rotazione che trasforma il piano in cui la frattura è contenuta nel piano xy, in modo che il primo lato sia parallelo all'asse x (è una scelta, se ne poteva fare anche un'altra)
inline Matrix<double,3,3> allinea_xy(const std::vector<Vector3d>& vertici, Vector3d normale)
{
    Matrix<double,3,3> rotazione;
    Vector3d v=(vertici[0]-vertici[1]).normalized();
    rotazione.row(0)=v;
    rotazione.row(1)=normale.cross(v);
    rotazione.row(2)=normale;
    return rotazione;
}




//funzione effettiva che memorizza le tracce nella struttura DFN
void memorizza_tracce(DFN& disc_frac_net,double tol)
{
    //controllo che la tolleranza scelta dall'utente non sia minore della precisione di macchina
    tol=std::max(tol,10*__DBL_EPSILON__);

    //--------- passo 1 ---------
    //scrematura grossolana delle fratture che sicuramente non si intersecano
    disc_frac_net.versori.reserve(disc_frac_net.numFratture);
    for(auto& poligono:disc_frac_net.vertici)
    {
        Vector3d norm=versore_normale(poligono);
        disc_frac_net.versori.push_back(norm);
    }
    std::vector<std::array<unsigned int,2>> fratture_suscettibili=scarta_fratture(disc_frac_net,tol);

    //--------- passo 2 ---------
    //riservo spazio per memorizzare i dati sulle varie strutture dati
    disc_frac_net.idTracce.reserve(fratture_suscettibili.size());
    disc_frac_net.estremiTracce.reserve(fratture_suscettibili.size());
    disc_frac_net.tips.reserve(fratture_suscettibili.size());
    disc_frac_net.tracce.reserve(fratture_suscettibili.size());
    disc_frac_net.lunghezze.reserve(fratture_suscettibili.size());
    unsigned int idt=0; //id delle tracce trovate

    //--------- passo 3 ---------
    //ricerca degli estremi della traccia
    for(auto& coppia:fratture_suscettibili)
    {
        //--------- output di ogni intersezione ---------
        std::array<unsigned int,2> tip={0,0};
        std::array<Vector3d,2> est_tra={Vector3d(NAN,NAN,NAN),Vector3d(NAN,NAN,NAN)};
        unsigned int pos=0; //è l'indice di riempimento di est_tra, ossia indica la posizione libera di est_tra (se pos==2 vuol dire che è pieno)

        //--------- passo 4 ---------
        //ruoto i poligoni sul piano xy per rendere più facili i confronti dopo
        Matrix<double,3,3> mat_rot=allinea_xy(disc_frac_net.vertici[coppia[0]],disc_frac_net.versori[coppia[0]]);
        std::vector<Vector3d> vertici_ruotati_2={}; //seconda frattura
        vertici_ruotati_2.resize(disc_frac_net.numVertici[coppia[1]]);

        for(unsigned int i=0;i<disc_frac_net.numVertici[coppia[1]];i++)
        {
            vertici_ruotati_2[i]=mat_rot*(disc_frac_net.vertici[coppia[1]][i]-disc_frac_net.vertici[coppia[0]][0]);
        }

        //--------- passo 5 ---------
        //triangolo la prima frattura della coppia
        std::vector<std::array<unsigned int,3>> triangolazione=triangola_frattura(disc_frac_net,coppia[0]);

        //--------- passo 6 ---------
        //ciclo su tutti i lati della seconda frattura della coppia per vedere se intersecano la prima frattura
        for(unsigned int i=0;i<disc_frac_net.numVertici[coppia[1]];i++)
        {
            //std::cout << i <<std::endl;
            if(!(vertici_ruotati_2[i](2)*vertici_ruotati_2[(i+1)%disc_frac_net.numVertici[coppia[1]]](2)>tol)) //controlla che il lato intersechi effettivamente il piano che contiene la frattura 1
            {
                //--------- passo 7 ---------
                //ciclo su tutti i triangoli della prima frattura nella coppia
                for(auto& triangolo:triangolazione)
                {
                    Matrix<double,3,3> mat_coeff;
                    Vector3d b=disc_frac_net.vertici[coppia[1]][i]-disc_frac_net.vertici[coppia[0]][triangolo[0]]; //vettore dei termini noti
                    mat_coeff.col(0)=disc_frac_net.vertici[coppia[0]][triangolo[1]]-disc_frac_net.vertici[coppia[0]][triangolo[0]];
                    mat_coeff.col(1)=disc_frac_net.vertici[coppia[0]][triangolo[2]]-disc_frac_net.vertici[coppia[0]][triangolo[0]];
                    mat_coeff.col(2)=disc_frac_net.vertici[coppia[1]][i]-disc_frac_net.vertici[coppia[1]][(i+1)%disc_frac_net.numVertici[coppia[1]]];
                    Vector3d x=mat_coeff.lu().solve(b);

                    //--------- passo 8 ---------
                    //controllo se la soluzione trovata appartiene al triangolo
                    if(x(0)>=-tol && x(1)>=-tol && x(0)+x(1)<=1.+tol && x(2)<=1.+tol && x(2)>=-tol) //l'ultimo controllo non dovrebbe servire ma per ora lo tengo
                    {
                        Vector3d y=0.5*(disc_frac_net.vertici[coppia[1]][i]+x(2)*(-mat_coeff.col(2))+disc_frac_net.vertici[coppia[0]][triangolo[0]]+x(0)*mat_coeff.col(0)+x(1)*mat_coeff.col(1));
                        est_tra[pos]=y;
                        tip[1]+=1;

                        //--------- passo 9 ---------
                        //controllo che il nuovo estremo non sia uguale all'eventuale estremo già trovato precedente
                        if(pos!=0)//devo fare 2 in innestati in quanto c'è comportamento indefinito quando pos=0 e faccio est_tra[pos-1]
                        {
                            if(!((y-est_tra[pos-1]).norm()/std::max(std::max(y.norm(),est_tra[pos-1].norm()),1.)>tol))
                            {
                                est_tra[pos]=Vector3d(NAN,NAN,NAN);
                                pos+=-1;
                                tip[1]+=-1;
                            }
                        }
                        //--------- passo 10 ---------
                        //controllo se la soluzione trovata appartiene al bordo del primo poligono nella coppia
                        if((x(0)+x(1)>=1.-tol && x(0)+x(1)<=1.+tol) || (x(0)<=1.+tol && !(abs(x(1))>tol) && (triangolo[1]==1)) || (x(1)<=1.+tol && !(abs(x(0))>tol) && (triangolo[2]==disc_frac_net.numVertici[coppia[0]]-1)))
                        {
                            tip[0]+=1;
                        }
                        pos++;
                        break; //dato che ho trovato che un estremo della traccia appartiene ad un triangolo, non devo controllare che appartenga ai triangoli successivi
                    }
                }
            }
            if(pos==2) //se ho trovato tutti e 2 gli estremi della traccia, ho finito
            {
                break;
            }
        }
        //--------- passo 11 ---------
        //se non ho trovato tutti e 2 gli estremi della traccia, vado avanti e ciclo sui lati della prima frattura della coppia e sui triangoli della seconda frattura della coppia per trovare l'altro/gli estremo/i della traccia
        if(pos!=2)
        {
            //non ricommento tutto in quanto è molto simile a prima (però non è proprio uguale identico, ci sono alcuni passaggi che qua non servono più)
            mat_rot=allinea_xy(disc_frac_net.vertici[coppia[1]],disc_frac_net.versori[coppia[0]]);
            std::vector<Vector3d> vertici_ruotati_1;
            vertici_ruotati_1.resize(disc_frac_net.numVertici[coppia[0]]);
            for(unsigned int i=0;i<disc_frac_net.numVertici[coppia[0]];i++)
            {
                vertici_ruotati_1[i]=mat_rot*(disc_frac_net.vertici[coppia[0]][i]-disc_frac_net.vertici[coppia[1]][0]);
            }

            triangolazione=triangola_frattura(disc_frac_net,coppia[1]);

            for(unsigned int i=0;i<disc_frac_net.numVertici[coppia[0]];i++)
            {
                if(!(vertici_ruotati_1[i](2)*vertici_ruotati_1[(i+1)%disc_frac_net.numVertici[coppia[0]]](2)>tol))
                {
                    for(auto& triangolo:triangolazione)
                    {
                        Matrix<double,3,3> mat_coeff;
                        Vector3d b=disc_frac_net.vertici[coppia[0]][i]-disc_frac_net.vertici[coppia[1]][triangolo[0]];
                        mat_coeff.col(0)=disc_frac_net.vertici[coppia[1]][triangolo[1]]-disc_frac_net.vertici[coppia[1]][triangolo[0]];
                        mat_coeff.col(1)=disc_frac_net.vertici[coppia[1]][triangolo[2]]-disc_frac_net.vertici[coppia[1]][triangolo[0]];
                        mat_coeff.col(2)=disc_frac_net.vertici[coppia[0]][i]-disc_frac_net.vertici[coppia[0]][(i+1)%disc_frac_net.numVertici[coppia[0]]];
                        Vector3d x=mat_coeff.lu().solve(b);
                        if(x(0)>=-tol && x(1)>=-tol && x(0)+x(1)<=1.+tol && x(2)<=1.+tol && x(2)>=-tol)
                        {
                            Vector3d y=0.5*(disc_frac_net.vertici[coppia[0]][i]+x(2)*(-mat_coeff.col(2))+disc_frac_net.vertici[coppia[1]][triangolo[0]]+x(0)*mat_coeff.col(0)+x(1)*mat_coeff.col(1));
                            est_tra[pos]=y;
                            tip[0]+=1;
                            if(pos!=0)
                            {
                                if(!((y-est_tra[pos-1]).norm()/std::max(std::max(y.norm(),est_tra[pos-1].norm()),1.)>tol))
                                {
                                    est_tra[pos]=Vector3d(NAN,NAN,NAN);
                                    pos+=-1;
                                    tip[0]+=-1;
                                }
                            }
                            if((x(0)+x(1)>=1.-tol && x(0)+x(1)<=1.+tol) || (x(0)<=1.+tol && !(abs(x(1))>tol) && (triangolo[1]==1)) || (x(1)<=1.+tol && !(abs(x(0))>tol) && (triangolo[2]==disc_frac_net.numVertici[coppia[1]]-1)))
                            {
                                tip[0]+=1;
                            }
                            pos++;
                            break;
                        }
                    }
                }
                if(pos==2)
                {
                    break;
                }
            }
        }
        //--------- passo 12 ---------
        //memorizzo i risultati e stampo un output dei risultati trovati
        if(pos==2) //se effettivamente ho trovato tutti e 2 gli estremi della traccia e non solo un punto di intersezione, che non contiamo
        {
            disc_frac_net.idTracce.push_back(idt);
            std::cout << "TRACCIA "<<idt<<":"<<std::endl;
            std::array<unsigned int,2> coppia2={disc_frac_net.idFratture[coppia[0]],disc_frac_net.idFratture[coppia[1]]};
            disc_frac_net.tracce.push_back(coppia2);
            std::cout << "la traccia e' stata generata dalle fratture "<<coppia2[0]<<" e "<<coppia2[1]<<std::endl;
            disc_frac_net.estremiTracce.push_back(est_tra);
            std::cout << "gli estremi della traccia sono "<<std::endl<<est_tra[0] <<std::endl<<"e"<<std::endl<<est_tra[1]<<std::endl;
            std::array<bool,2> tip2={(tip[0]==2),(tip[1]==2)};
            disc_frac_net.tips.push_back(tip2);
            std::cout << "la traccia e' passante per la frattura "<<coppia[0]<<"? "<<tip2[0]<<" (0=falso, 1=vero)"<<std::endl;
            std::cout << "la traccia e' passante per la frattura "<<coppia[1]<<"? "<<tip2[1]<<" (0=falso, 1=vero)"<<std::endl;
            disc_frac_net.lunghezze.push_back((est_tra[0]-est_tra[1]).norm());
            std::cout << "la lunghezza della traccia e' "<<(est_tra[0]-est_tra[1]).norm()<<std::endl;
            std::cout << "-----------------------------------------------------"<<std::endl;
            idt++;
        }
    }
    disc_frac_net.numTracce=idt;
    std::cout << "Alla fine sono state trovate "<<idt<<" tracce"<<std::endl;
}




//***********************************************PARTE 2***********************************************




//funzione che trova l'intersezione tra 2 segmenti, in verità trova il punto medio del segmento di lunghezza minima che collega i 2 segmenti, che, se i 2 segmenti sono complanari, coincide con il punto di intersezione
std::tuple<Vector3d,double> interseca_segmenti(const Vector3d& A1,const Vector3d& A2,const Vector3d& B1,const Vector3d& B2)
{
    Vector3d n=(A2-A1).normalized();
    Vector3d v=B2-B1;
    Vector3d w=B1-A1;
    double alfa=(v.dot(w.dot(n)*n-w))/(v.squaredNorm()-(v.dot(n))*(v.dot(n)));

    n=(B2-B1).normalized();
    v=A2-A1;
    w=-w;
    double beta=(v.dot(w.dot(n)*n-w))/(v.squaredNorm()-(v.dot(n))*(v.dot(n)));
    Vector3d x=0.5*(B1+alfa*(B2-B1)+A1+beta*(A2-A1));
    return std::make_tuple(x,beta);
}




//funzione che divide il poligono in 2 sottopoligoni in base al segmento tagliante e di conseguenza aggiorna la mesh
void aggiorna_mesh(PolygonalMesh& mesh
                   ,const std::vector<unsigned int>& poligoni
                   ,const std::vector<std::array<unsigned int,2>>& lati_coinvolti
                   ,const Vector3d& est1
                   ,const Vector3d& est2
                   ,const double& tol)
{
    //------PASSO 1------: per ogni poligono coinvolto nel taglio aggiorno le informazioni relative alle celle 0D e 1D già esistenti, considerando il nuovo taglio
    std::unordered_map<unsigned int,std::vector<unsigned int>> mappa_vecchi_lati_nuovi_lati; // chiave: lato vecchio | contenuto: coppia di lati che sostituiscono il lato vecchio
    for(unsigned int i=0; i<poligoni.size();i++)
    {
        std::array<unsigned int,2> est_taglio; //id estremi taglio
        unsigned int pos=0;
        //NOTA: per costruzione, se pos==2, non succederà che il programma tenterà di scrivere a est_taglio[pos], generando un comportamento indefinito

        //------PASSO 2------: controllo che la traccia non sia contenuta in una cella 2D
        if(lati_coinvolti[i][1]==std::numeric_limits<unsigned int>::max())
        {
            continue; //la traccia non mi va a modificare la mesh
        }

        //------PASSO 3------: aggiorno i punti d'intersezione a meno che questi non esistano già
        for(const unsigned int& lato:lati_coinvolti[i])
        {
            if(mappa_vecchi_lati_nuovi_lati.find(lato)==mappa_vecchi_lati_nuovi_lati.end()) //se l'intersezione con il lato non è stata ancora inserita
            {
                Vector3d& A1=mesh.Cell0DCoordinates[mesh.Cell1DVertices[lato][0]];
                Vector3d& A2=mesh.Cell0DCoordinates[mesh.Cell1DVertices[lato][1]];
                std::tuple<Vector3d,double> x=interseca_segmenti(A1,A2,est1,est2);
                if(std::abs(std::get<double>(x))>tol) //se il vertice di intersezione non è ugugale al primo estremo del lato
                {
                    if(std::abs(std::get<double>(x)-1)>tol) //se il vertice di intersezione non è uguale al secondo estremo del lato
                    {
                        //aggiorno le celle 0D
                        unsigned int ID_0D=mesh.NumberCell0D;
                        mesh.NumberCell0D++;
                        mesh.Cell0DId.push_back(ID_0D);
                        mesh.Cell0DCoordinates.push_back(std::get<Vector3d>(x));
                        mesh.Cell0DMarkers.push_back(mesh.Cell1DMarkers[lato]); //il marker del nuovo vertice viene ereditato dal marker del lato a cui appartiene

                        //aggiorno le celle 1D
                        unsigned int ID_1D=mesh.NumberCell1D;
                        mesh.NumberCell1D++;
                        mesh.Cell1DId.push_back(ID_1D);
                        std::array<unsigned int,2> vertici={ID_0D,mesh.Cell1DVertices[lato][1]};
                        mesh.Cell1DVertices[lato]={ID_0D,mesh.Cell1DVertices[lato][0]};
                        mesh.Cell1DVertices.push_back(vertici);
                        mesh.Cell1DMarkers.push_back(mesh.Cell1DMarkers[lato]);

                        //metto la nuova info nel dizionario e aggiorno est_taglio
                        mappa_vecchi_lati_nuovi_lati[lato]={lato,ID_1D};
                        est_taglio[pos]=ID_0D;
                    }
                    else //il vertice d'intersezione è uguale al secondo estremo del lato
                    {
                        est_taglio[pos]=mesh.Cell1DVertices[lato][1];
                        mappa_vecchi_lati_nuovi_lati[lato]={lato};
                    }
                }
                else //il vertice d'intersezione è uguale al primo estremo del lato
                {
                    est_taglio[pos]=mesh.Cell1DVertices[lato][0];
                    mappa_vecchi_lati_nuovi_lati[lato]={lato};
                }
            }
            else //l'intersezione dell'asse che contiene la traccia con il lato esiste già
            {
                // da come sono stati implementati i dati nella mappa, questo è un modo per recuperare l'ID dell'estremo della traccia che giace sul lato
                est_taglio[pos]=mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[lato][0]][0];
            }
            pos++;
        }

        //------PASSO 4------: aggiorno le altre celle 2d che hanno lati in comune con il poligono coinvolto considerato un lato coinvolto nel taglio
        if(mappa_vecchi_lati_nuovi_lati[lati_coinvolti[i][0]].size()==2 && mappa_vecchi_lati_nuovi_lati[lati_coinvolti[i][1]].size()==2) //faccio i 2 controlli in una volta sola
        {
            for(auto& cella2D:mesh.Cell2DId) //controllo su tutte le celle 2D
            {
                if(std::find(poligoni.begin(),poligoni.end(),cella2D)==poligoni.end() && cella2D<mesh.NumberCell2D-poligoni.size()) //devo controllare sui poligoni che non sono coinvolti nel taglio
                {
                    auto trova_vecchio_lato_coinvolto_0=std::find(mesh.Cell2DEdges[cella2D].begin(),mesh.Cell2DEdges[cella2D].end(),lati_coinvolti[i][0]);
                    unsigned int posizione_vecchio_lato_coinvolto_0=std::distance(mesh.Cell2DEdges[cella2D].begin(),trova_vecchio_lato_coinvolto_0);
                    if(trova_vecchio_lato_coinvolto_0!=mesh.Cell2DEdges[cella2D].end()) //se trovo il vecchio_lato_coinvolto_0 come lato del poligono
                    {
                        //aggiorno i vertici
                        unsigned int successivo=(posizione_vecchio_lato_coinvolto_0+1)%(mesh.Cell2DEdges[cella2D].size());
                        mesh.Cell2DVertices[cella2D].insert(mesh.Cell2DVertices[cella2D].begin()+successivo,mesh.Cell1DVertices[lati_coinvolti[i][0]][0]);
                        //aggiorno i lati
                        if(mesh.Cell1DVertices[lati_coinvolti[i][0]][1]==mesh.Cell2DVertices[cella2D][posizione_vecchio_lato_coinvolto_0])
                        {
                            mesh.Cell2DEdges[cella2D].insert(mesh.Cell2DEdges[cella2D].begin()+successivo,mappa_vecchi_lati_nuovi_lati[lati_coinvolti[i][0]][1]);
                        }
                        else
                        {
                            auto& precedente=trova_vecchio_lato_coinvolto_0;
                            mesh.Cell2DEdges[cella2D].insert(precedente,mappa_vecchi_lati_nuovi_lati[lati_coinvolti[i][0]][1]);
                        }
                    }
                    auto trova_vecchio_lato_coinvolto_1=std::find(mesh.Cell2DEdges[cella2D].begin(),mesh.Cell2DEdges[cella2D].end(),lati_coinvolti[i][1]);
                    unsigned int posizione_vecchio_lato_coinvolto_1=std::distance(mesh.Cell2DEdges[cella2D].begin(),trova_vecchio_lato_coinvolto_1);
                    if(trova_vecchio_lato_coinvolto_1!=mesh.Cell2DEdges[cella2D].end())
                    {
                        //aggiorno i vertici
                        unsigned int successivo=(posizione_vecchio_lato_coinvolto_1+1)%(mesh.Cell2DEdges[cella2D].size());
                        mesh.Cell2DVertices[cella2D].insert(mesh.Cell2DVertices[cella2D].begin()+successivo,mesh.Cell1DVertices[lati_coinvolti[i][1]][0]);
                        //aggiorno i lati
                        if(mesh.Cell1DVertices[lati_coinvolti[i][1]][1]==mesh.Cell2DVertices[cella2D][posizione_vecchio_lato_coinvolto_1])
                        {
                            mesh.Cell2DEdges[cella2D].insert(mesh.Cell2DEdges[cella2D].begin()+successivo,mappa_vecchi_lati_nuovi_lati[lati_coinvolti[i][1]][1]);
                        }
                        else
                        {
                            auto& precedente=trova_vecchio_lato_coinvolto_1;
                            mesh.Cell2DEdges[cella2D].insert(precedente,mappa_vecchi_lati_nuovi_lati[lati_coinvolti[i][1]][1]);
                        }
                    }
                }
            }
        }
        if(mappa_vecchi_lati_nuovi_lati[lati_coinvolti[i][0]].size()==2 && mappa_vecchi_lati_nuovi_lati[lati_coinvolti[i][1]].size()==1) //faccio 1 controllo solo
        {
            for(auto& cella2D:mesh.Cell2DId) //controllo su tutte le celle 2D
            {
                if(std::find(poligoni.begin(),poligoni.end(),cella2D)==poligoni.end() && cella2D<mesh.NumberCell2D-poligoni.size()) //devo controllare sui poligoni che non sono coinvolti nel taglio
                {
                    auto trova_vecchio_lato_coinvolto_0=std::find(mesh.Cell2DEdges[cella2D].begin(),mesh.Cell2DEdges[cella2D].end(),lati_coinvolti[i][0]);
                    unsigned int posizione_vecchio_lato_coinvolto_0=std::distance(mesh.Cell2DEdges[cella2D].begin(),trova_vecchio_lato_coinvolto_0);
                    if(trova_vecchio_lato_coinvolto_0!=mesh.Cell2DEdges[cella2D].end()) //se trovo il vecchio_lato_coinvolto_0 come lato del poligono
                    {
                        //aggiorno i vertici
                        unsigned int successivo=(posizione_vecchio_lato_coinvolto_0+1)%(mesh.Cell2DEdges[cella2D].size());
                        mesh.Cell2DVertices[cella2D].insert(mesh.Cell2DVertices[cella2D].begin()+successivo,mesh.Cell1DVertices[lati_coinvolti[i][0]][0]);
                        //aggiorno i lati
                        if(mesh.Cell1DVertices[lati_coinvolti[i][0]][1]==mesh.Cell2DVertices[cella2D][posizione_vecchio_lato_coinvolto_0])
                        {
                            mesh.Cell2DEdges[cella2D].insert(mesh.Cell2DEdges[cella2D].begin()+successivo,mappa_vecchi_lati_nuovi_lati[lati_coinvolti[i][0]][1]);
                        }
                        else
                        {
                            auto& precedente=trova_vecchio_lato_coinvolto_0;
                            mesh.Cell2DEdges[cella2D].insert(precedente,mappa_vecchi_lati_nuovi_lati[lati_coinvolti[i][0]][1]);
                        }
                    }
                }
            }
        }
        if(mappa_vecchi_lati_nuovi_lati[lati_coinvolti[i][0]].size()==1 && mappa_vecchi_lati_nuovi_lati[lati_coinvolti[i][1]].size()==2) //faccio 1 controllo solo
        {
            for(auto& cella2D:mesh.Cell2DId) //controllo su tutte le celle 2D
            {
                if(std::find(poligoni.begin(),poligoni.end(),cella2D)==poligoni.end() && cella2D<mesh.NumberCell2D-poligoni.size()) //devo controllare sui poligoni che non sono coinvolti nel taglio
                {
                    auto trova_vecchio_lato_coinvolto_1=std::find(mesh.Cell2DEdges[cella2D].begin(),mesh.Cell2DEdges[cella2D].end(),lati_coinvolti[i][1]);
                    unsigned int posizione_vecchio_lato_coinvolto_1=std::distance(mesh.Cell2DEdges[cella2D].begin(),trova_vecchio_lato_coinvolto_1);
                    if(trova_vecchio_lato_coinvolto_1!=mesh.Cell2DEdges[cella2D].end())
                    {
                        //aggiorno i vertici
                        unsigned int successivo=(posizione_vecchio_lato_coinvolto_1+1)%(mesh.Cell2DEdges[cella2D].size());
                        mesh.Cell2DVertices[cella2D].insert(mesh.Cell2DVertices[cella2D].begin()+successivo,mesh.Cell1DVertices[lati_coinvolti[i][1]][0]);
                        //aggiorno i lati
                        if(mesh.Cell1DVertices[lati_coinvolti[i][1]][1]==mesh.Cell2DVertices[cella2D][posizione_vecchio_lato_coinvolto_1])
                        {
                            mesh.Cell2DEdges[cella2D].insert(mesh.Cell2DEdges[cella2D].begin()+successivo,mappa_vecchi_lati_nuovi_lati[lati_coinvolti[i][1]][1]);
                        }
                        else
                        {
                            auto& precedente=trova_vecchio_lato_coinvolto_1;
                            mesh.Cell2DEdges[cella2D].insert(precedente,mappa_vecchi_lati_nuovi_lati[lati_coinvolti[i][1]][1]);
                        }
                    }
                }
            }
        }

        //------PASSO 5------: aggiungo le info relative al taglio
        //celle 1D
        unsigned int ID_1D=mesh.NumberCell1D;
        mesh.Cell1DId.push_back(ID_1D);
        mesh.NumberCell1D++;
        mesh.Cell1DMarkers.push_back(0); //il marker sicuramente è nullo dato che il segmento taglia la frattura (il caso degenere è stato già trattato prima)
        mesh.Cell1DVertices.push_back(est_taglio);

        //celle 2D
        mesh.Cell2DId.push_back(mesh.NumberCell2D);
        mesh.NumberCell2D++;

        //------PASSO 6------: creazione del primo poligono
        std::array<std::vector<unsigned int>,2> poligono=nuovo_poligono(mesh,mappa_vecchi_lati_nuovi_lati,poligoni,i,lati_coinvolti[i][0],lati_coinvolti[i][1],est_taglio[0],est_taglio[1],ID_1D);
        std::vector<unsigned int> errore={std::numeric_limits<unsigned int>::max()};
        if(poligono[0]!=errore)
        {
            mesh.Cell2DEdges.push_back(poligono[1]);
            mesh.Cell2DVertices.push_back(poligono[0]);
        }
        else
        {
            std::cerr << "Non e' stato possibile aggiornare i nuovi poligoni dato che c'e' un errore"<<std::endl;
            return;
        }
        //------PASSO 7------: creazione del secondo poligono
        poligono=nuovo_poligono(mesh,mappa_vecchi_lati_nuovi_lati,poligoni,i,lati_coinvolti[i][1],lati_coinvolti[i][0],est_taglio[1],est_taglio[0],ID_1D);
        if(poligono[0]!=errore)
        {
            mesh.Cell2DEdges[poligoni[i]]=poligono[1];
            mesh.Cell2DVertices[poligoni[i]]=poligono[0];
        }
        else
        {
            std::cerr << "Non e' stato possibile aggiornare i nuovi poligoni dato che c'e' un errore"<<std::endl;
            return;
        }

    }
}



std::array<std::vector<unsigned int>,2> nuovo_poligono(PolygonalMesh& mesh,
                                                       std::unordered_map<unsigned int, std::vector<unsigned int>>& mappa_vecchi_lati_nuovi_lati,
                                                       const std::vector<unsigned int>& poligoni,
                                                       const unsigned int& i,
                                                       const unsigned int& vecchio_lato_partenza,
                                                       const unsigned int& vecchio_lato_arrivo,
                                                       const unsigned int& est_taglio_partenza,
                                                       const unsigned int& est_taglio_arrivo,
                                                       const unsigned int& ID_1D)
{
    std::vector<unsigned int> vertici_nuovo_poligono;
    std::vector<unsigned int> lati_nuovo_poligono;
    unsigned int larghezza=mesh.Cell2DEdges[poligoni[i]].size();
    vertici_nuovo_poligono.reserve(larghezza);
    lati_nuovo_poligono.reserve(larghezza);
    //vado dal primo lato coinvolto al secondo lato coinvolto
    vertici_nuovo_poligono.push_back(est_taglio_arrivo);
    vertici_nuovo_poligono.push_back(est_taglio_partenza);
    lati_nuovo_poligono.push_back(ID_1D);
    auto trova_lato_vecchio_1=std::find(mesh.Cell2DEdges[poligoni[i]].begin(),mesh.Cell2DEdges[poligoni[i]].end(),vecchio_lato_partenza);
    unsigned int posto=(std::distance(mesh.Cell2DEdges[poligoni[i]].begin(),trova_lato_vecchio_1)+1)%larghezza;
    bool finito=false;
    bool mancano_info=false;
    //Vanno considerati vari casi
    if(mesh.Cell2DEdges[poligoni[i]][posto]!=vecchio_lato_arrivo) //caso normale
    {
        if(mesh.Cell2DVertices[poligoni[i]][posto]!=est_taglio_partenza)
        {
            vertici_nuovo_poligono.push_back(mesh.Cell2DVertices[poligoni[i]][posto]);
            lati_nuovo_poligono.push_back(0);
            mancano_info=true;
        }
        else
        {
            lati_nuovo_poligono.push_back(mesh.Cell2DEdges[poligoni[i]][posto]);
            posto++;
            posto=posto%larghezza;
        }
    }
    else //caso degenere: ho un triangolo, che va gestito in una maniera diversa
    {
        //Siano A,B i nuovi lati derivati da vecchio_lato_partenza e C,D i nuovi lati derivati da vecchio_lato_arrivo se esistono tutti quanti
        bool AcongC=(std::find(mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_arrivo][0]].begin(),mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_arrivo][0]].end(),mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_partenza][0]][0])!=mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_arrivo][0]].end())
                 || (std::find(mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_arrivo][0]].begin(),mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_arrivo][0]].end(),mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_partenza][0]][1])!=mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_arrivo][0]].end());
        if(AcongC) // A congiunge C?
        {
            vertici_nuovo_poligono.push_back(mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_partenza][0]][1]);
            lati_nuovo_poligono.push_back(mappa_vecchi_lati_nuovi_lati[vecchio_lato_partenza][0]);
            lati_nuovo_poligono.push_back(mappa_vecchi_lati_nuovi_lati[vecchio_lato_arrivo][0]);
        }
        else
        {
            if(mappa_vecchi_lati_nuovi_lati[vecchio_lato_arrivo].size()==2) //ho C e D?
            {
                bool AcongD=(std::find(mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_arrivo][1]].begin(),mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_arrivo][1]].end(),mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_partenza][0]][0])!=mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_arrivo][1]].end())
                         || (std::find(mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_arrivo][1]].begin(),mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_arrivo][1]].end(),mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_partenza][0]][1])!=mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_arrivo][1]].end());
                if(AcongD)
                {
                    vertici_nuovo_poligono.push_back(mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_partenza][0]][1]);
                    lati_nuovo_poligono.push_back(mappa_vecchi_lati_nuovi_lati[vecchio_lato_partenza][0]);
                    lati_nuovo_poligono.push_back(mappa_vecchi_lati_nuovi_lati[vecchio_lato_arrivo][1]);
                }
                else // B congiunge D
                {
                    if(mappa_vecchi_lati_nuovi_lati[vecchio_lato_arrivo].size()==1)
                    {
                        std::cerr<<"C'e' un errore nel taglio del poligono "<< poligoni[i]<<std::endl;
                        vertici_nuovo_poligono={std::numeric_limits<unsigned int>::max()};
                        lati_nuovo_poligono={std::numeric_limits<unsigned int>::max()};
                        return {vertici_nuovo_poligono,lati_nuovo_poligono};
                    }
                    vertici_nuovo_poligono.push_back(mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_partenza][1]][1]);
                    lati_nuovo_poligono.push_back(mappa_vecchi_lati_nuovi_lati[vecchio_lato_partenza][1]);
                    lati_nuovo_poligono.push_back(mappa_vecchi_lati_nuovi_lati[vecchio_lato_arrivo][1]);
                }
            }
            else // B congiunge C
            {
                if(mappa_vecchi_lati_nuovi_lati[vecchio_lato_arrivo].size()==1)
                {
                    std::cerr<<"C'e' un errore nel taglio del poligono "<< poligoni[i]<<std::endl;
                    vertici_nuovo_poligono={std::numeric_limits<unsigned int>::max()};
                    lati_nuovo_poligono={std::numeric_limits<unsigned int>::max()};
                    return {vertici_nuovo_poligono,lati_nuovo_poligono};
                }
                vertici_nuovo_poligono.push_back(mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_partenza][1]][1]);
                lati_nuovo_poligono.push_back(mappa_vecchi_lati_nuovi_lati[vecchio_lato_partenza][1]);
                lati_nuovo_poligono.push_back(mappa_vecchi_lati_nuovi_lati[vecchio_lato_arrivo][0]);
            }
        }
        finito=true;
    }
    if(!finito) //nel caso normale devo ancora finire
    {
        if(mancano_info)
        {
            std::array<unsigned int,2> lato_1={est_taglio_partenza,vertici_nuovo_poligono.back()};
            auto trova_lato_1=std::find(mesh.Cell1DVertices.begin(),mesh.Cell1DVertices.end(),lato_1);
            if(trova_lato_1==mesh.Cell1DVertices.end())
            {
                lato_1={vertici_nuovo_poligono.back(),est_taglio_partenza};
                trova_lato_1=std::find(mesh.Cell1DVertices.begin(),mesh.Cell1DVertices.end(),lato_1);
            }
            unsigned int posizione=std::distance(mesh.Cell1DVertices.begin(),trova_lato_1);
            lati_nuovo_poligono[1]=mesh.Cell1DId[posizione];
            lati_nuovo_poligono.push_back(mesh.Cell2DEdges[poligoni[i]][posto]);
            posto++;
            posto=posto%larghezza;
            mancano_info=false;
        }
        while(mesh.Cell2DEdges[poligoni[i]][posto]!=vecchio_lato_arrivo)
        {
            vertici_nuovo_poligono.push_back(mesh.Cell2DVertices[poligoni[i]][posto]);
            lati_nuovo_poligono.push_back(mesh.Cell2DEdges[poligoni[i]][posto]);
            posto++;
            posto=posto%larghezza;
        }
        if(mesh.Cell2DVertices[poligoni[i]][posto]!=est_taglio_arrivo)
        {
            vertici_nuovo_poligono.push_back(mesh.Cell2DVertices[poligoni[i]][posto]);
            lati_nuovo_poligono.push_back(0);
            mancano_info=true;
        }
        if(mancano_info)
        {
            std::array<unsigned int,2> lato_1={est_taglio_arrivo,vertici_nuovo_poligono.back()};
            auto trova_lato_1=std::find(mesh.Cell1DVertices.begin(),mesh.Cell1DVertices.end(),lato_1);
            if(trova_lato_1==mesh.Cell1DVertices.end())
            {
                lato_1={vertici_nuovo_poligono.back(),est_taglio_arrivo};
                trova_lato_1=std::find(mesh.Cell1DVertices.begin(),mesh.Cell1DVertices.end(),lato_1);
            }
            unsigned int posizione=std::distance(mesh.Cell1DVertices.begin(),trova_lato_1);
            lati_nuovo_poligono.back()=mesh.Cell1DId[posizione];
        }
    }
    return {vertici_nuovo_poligono,lati_nuovo_poligono};
}




//funzione che trova i poligoni che intersecano il segmento
//OUTPUT: tupla di 2 elementi: 1) vettore degli indici dei poligoni intersecati dal segmento 2)vettori di vettori dei lati che intersecano l'asse che contiene il segmento
std::tuple<std::vector<unsigned int>,std::vector<std::array<unsigned int,2>>> contatto_poligoni_segmento(const DFN& disc_frac_net,const unsigned int& id,PolygonalMesh& mesh,const Vector3d& est1,const Vector3d& est2, const double& tol)
{
    //riga 1 e riga 3 della matrice di rotazione per spostare il problema sul piano yz (la riga 2 non serve)
    Vector3d vec1_rot=est2-est1;
    double norma=vec1_rot.norm();
    vec1_rot=1./norma*vec1_rot;
    Vector3d vec2_rot=(vec1_rot.cross(disc_frac_net.versori[id]));

    std::vector<unsigned int> temp_pol; //output 1
    std::vector<std::array<unsigned int,2>> temp_lat; //output 2
    temp_pol.reserve(mesh.NumberCell2D);
    temp_lat.reserve(mesh.NumberCell2D);

    //funzione lambda che calcola l'intersezione tra il lato e l'asse che contiene la traccia
    auto ottieni_intersezione = [&](const Vector3d& p1, const Vector3d& p2) -> double
    {
        double x1=vec1_rot.dot(p1-est1);
        double z1=vec2_rot.dot(p1-est1);
        double x2=vec1_rot.dot(p2-est1);
        double z2=vec2_rot.dot(p2-est1);
        if (!(z2*z1>tol))
        {
            return x2-z2*(x1-x2)/(z1-z2);
        }
        return NAN;
    };

    //selezione dei poligoni che intersecano il segmento
    for(unsigned int n2d=0; n2d < mesh.NumberCell2D; n2d++) //ciclo su tutti i poligoni della mesh
    {
        //vettore dei lati che intersecano l'asse del segmento
        std::array<unsigned int,2> lati_coinvolti={std::numeric_limits<unsigned int>::max(),std::numeric_limits<unsigned int>::max()}; //è improbabile che ci sia veramente un lato coinvolto che abbia questo indice
        unsigned int pos=0;

        //estremi della sezione
        double min_sez=norma;
        double max_sez=0.;

        //serve dopo
        double intersezione_precedente=NAN;
        for(auto& lato:mesh.Cell2DEdges[n2d])
        {
            Vector3d& p1=mesh.Cell0DCoordinates[mesh.Cell1DVertices[lato][0]];
            Vector3d& p2=mesh.Cell0DCoordinates[mesh.Cell1DVertices[lato][1]];
            double intersezione=ottieni_intersezione(p1,p2);
            if(!std::isnan(intersezione)) //se c'è intersezione tra il lato e l'asse che contiene la traccia
            {
                if(!std::isnan(intersezione_precedente)) //se anche il lato prima si intersecava
                {
                    if(std::abs(intersezione-intersezione_precedente)/std::max(std::max(intersezione,intersezione_precedente),1.)>tol && pos!=2) //sto controllo serve quando ho l'intersezione della retta con un lato che avviene proprio in un estremo del lato e va così a coinvolgere anche il lato successivo
                    {
                        lati_coinvolti[pos]=lato;
                        pos++;
                        min_sez=std::min(std::max(std::min(min_sez,intersezione),0.),norma);
                        max_sez=std::min(std::max(std::max(max_sez,intersezione),0.),norma);
                    }
                }
                else //se invece il lato prima non s'intersecava, non serve fare il controllo
                {
                    lati_coinvolti[pos]=lato;
                    pos++;
                    min_sez=std::min(std::max(std::min(min_sez,intersezione),0.),norma);
                    max_sez=std::min(std::max(std::max(max_sez,intersezione),0.),norma);
                }
            }
            //aggiorno
            intersezione_precedente=intersezione;
        }
        //se il poligono interseca il segmento in un insieme di misura non nulla, salvo i dati
        if(max_sez-min_sez>tol)
        {
            temp_lat.push_back(lati_coinvolti);
            temp_pol.push_back(n2d);
        }
    }
    return std::make_tuple(temp_pol,temp_lat);
}




//procedura che crea le mesh sulle fratture a partire dalle tracce
void definisci_mesh(DFN& disc_frac_net, double tol)
{
    //controllo che la tolleranza scelta dall'utente non sia minore della precisione di macchina
    tol=std::max(tol,10*__DBL_EPSILON__);

    disc_frac_net.meshPoligonali.reserve(disc_frac_net.numFratture);
    for(unsigned int id=0;id<disc_frac_net.numFratture;id++) // per ogni frattura
    {
        PolygonalMesh mesh;

        //riservo memoria per gli attributi della struttura
        unsigned int num_tracce_per_frattura=disc_frac_net.traccePassanti[id].size()+disc_frac_net.tracceNonPassanti[id].size();
        unsigned int capacita=disc_frac_net.numVertici[id]+(num_tracce_per_frattura+2)*num_tracce_per_frattura; //nel caso peggiore ogni frattura si interseca con tutte le altre
        mesh.Cell0DId.reserve(capacita);
        mesh.Cell0DCoordinates.reserve(capacita);
        mesh.Cell0DMarkers.reserve(capacita);

        capacita=disc_frac_net.numVertici[id]+(num_tracce_per_frattura+1)*num_tracce_per_frattura; //nel caso peggiore ogni frattura si interseca con tutte le altre
        mesh.Cell1DId.reserve(capacita);
        mesh.Cell1DVertices.reserve(capacita);
        mesh.Cell1DMarkers.reserve(capacita);

        capacita=num_tracce_per_frattura*num_tracce_per_frattura; //nel caso peggiore ogni frattura si interseca con tutte le altre
        mesh.Cell2DId.reserve(capacita);
        mesh.Cell2DVertices.reserve(capacita);
        mesh.Cell2DEdges.reserve(capacita);

        //questo mi serve dopo
        std::vector<unsigned int> pol_vertici_lati;
        pol_vertici_lati.reserve(disc_frac_net.numVertici[id]);

        mesh.NumberCell0D=disc_frac_net.numVertici[id];
        mesh.NumberCell1D=disc_frac_net.numVertici[id];
        for(unsigned int i=0;i<disc_frac_net.numVertici[id];i++)
        {
            //celle0D
            mesh.Cell0DId.push_back(i);
            mesh.Cell0DCoordinates.push_back(disc_frac_net.vertici[id][i]);
            mesh.Cell0DMarkers.push_back(i);

            //celle1D
            mesh.Cell1DId.push_back(i);
            std::array<unsigned int,2> est_cella1D={i,(i+1)%disc_frac_net.numVertici[id]};
            mesh.Cell1DVertices.push_back(est_cella1D);
            mesh.Cell1DMarkers.push_back(i);

            pol_vertici_lati.push_back(i+disc_frac_net.numVertici[id]);
        }
        //celle2D
        mesh.Cell2DVertices.push_back(pol_vertici_lati);
        mesh.Cell2DEdges.push_back(pol_vertici_lati);
        mesh.Cell2DId={0};
        mesh.NumberCell2D=1;

        for(auto& traPas:disc_frac_net.traccePassanti[id]) //ciclo su tutte le tracce passanti
        {
            Vector3d est1=disc_frac_net.estremiTracce[traPas][0];
            Vector3d est2=disc_frac_net.estremiTracce[traPas][1];
            std::tuple<std::vector<unsigned int>,std::vector<std::array<unsigned int,2>>> poligoni_e_lati_coinvolti=contatto_poligoni_segmento(disc_frac_net,id,mesh,est1,est2,tol);
            aggiorna_mesh(mesh,std::get<std::vector<unsigned int>>(poligoni_e_lati_coinvolti),std::get<std::vector<std::array<unsigned int,2>>>(poligoni_e_lati_coinvolti),est1,est2,tol);
        }
        for(auto& traNoPas:disc_frac_net.tracceNonPassanti[id]) //ciclo su tutte le tracce non passanti
        {
            Vector3d est1=disc_frac_net.estremiTracce[traNoPas][0];
            Vector3d est2=disc_frac_net.estremiTracce[traNoPas][1];
            std::tuple<std::vector<unsigned int>,std::vector<std::array<unsigned int,2>>> poligoni_e_lati_coinvolti=contatto_poligoni_segmento(disc_frac_net,id,mesh,est1,est2,tol);
            aggiorna_mesh(mesh,std::get<std::vector<unsigned int>>(poligoni_e_lati_coinvolti),std::get<std::vector<std::array<unsigned int,2>>>(poligoni_e_lati_coinvolti),est1,est2,tol);
        }
        disc_frac_net.meshPoligonali.push_back(mesh);
    }
}
}
