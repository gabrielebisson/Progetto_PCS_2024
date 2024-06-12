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




//funzione che divide il poligono in 2 sottopoligoni in base al segmento tagliante
//assumo sempre che lati adiacenti non siano collineari
void aggiorna_mesh(PolygonalMesh& mesh
                   ,const std::vector<unsigned int>& poligoni
                   ,const std::vector<std::array<unsigned int,2>>& lati_coinvolti
                   ,const Vector3d& est1
                   ,const Vector3d& est2
                   ,const double& tol)
{
    std::map<unsigned int,std::vector<unsigned int>> mappa_vecchi_lati_nuovi_lati; //dizionario che associa ad ogni vecchio lato gli ID dei 2 nuovi lati che lo sostituiranno
    for(unsigned int i=0; i<poligoni.size();i++) //indago le modifiche su tutti i poligoni coinvolti dal taglio
    {
        //qui vado a salvare gli ID degli estremi del taglio rispettivi al poligono considerato
        std::array<unsigned int,2> est_taglio;
        unsigned int pos=0;
        //NOTA: per costruzione, se pos==2, non succederà che il programma tenterà di scrivere a est_taglio[pos], generando un comportamento indefinito

        //CONTROLLO PRELIMINARE: il segmento è contenuto in un lato del bordo?
        //NOTA: le tracce sono generate solo da 2 fratture (assunzione di base), quindi se dovessi trovare un array di lati coinvolti che come secondo membro ha NAN, vuol dire che la traccia è contenuta in un lato che sicuramente è di bordo
        if(lati_coinvolti[i][1]==std::numeric_limits<unsigned int>::max())
        {
            break; //non faccio assolutamente niente dato che la traccia non mi va a modificare la mesh
        }

        //vado ad aggiornare i punti d'intersezione a meno che questi non esistano già
        for(const unsigned int& lato:lati_coinvolti[i])
        {
            if(mappa_vecchi_lati_nuovi_lati.find(lato)==mappa_vecchi_lati_nuovi_lati.end()) //se l'intersezione con il lato non è stata ancora inserita
            {
                Vector3d& A1=mesh.Cell0DCoordinates[mesh.Cell1DVertices[lato][0]];
                Vector3d& A2=mesh.Cell0DCoordinates[mesh.Cell1DVertices[lato][1]];
                std::tuple<Vector3d,double> x=interseca_segmenti(A1,A2,est1,est2);
                if(std::abs(std::get<double>(x))>tol) //se il vertice di intersezione non è ugugale al primo estremo del lato
                {
                    if(std::abs(std::get<double>(x)-1)>tol) //se il vertice di intersezione non è ugugale al secondo estremo del lato
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
                        mesh.Cell1DVertices[lato]={mesh.Cell1DVertices[lato][0],ID_0D};
                        mesh.Cell1DVertices.push_back(vertici);
                        mesh.Cell1DMarkers.push_back(mesh.Cell1DMarkers[lato]);

                        //metto la nuova info nel dizionario e aggiorno est_taglio
                        mappa_vecchi_lati_nuovi_lati[lato]={lato,ID_1D};
                        est_taglio[pos]=ID_0D;
                        pos++;
                    }
                    else
                    {
                        unsigned int v1=mesh.Cell1DVertices[lato][0];
                        unsigned int v2=mesh.Cell1DVertices[lato][1];
                        mesh.Cell1DVertices[lato]={v2,v1};
                        est_taglio[pos]=mesh.Cell1DVertices[lato][0];
                        pos++;
                        mappa_vecchi_lati_nuovi_lati[lato]={lato};
                    }
                }
                else
                {
                    est_taglio[pos]=mesh.Cell1DVertices[lato][0];
                    pos++;
                    mappa_vecchi_lati_nuovi_lati[lato]={lato};
                }
            }
            else
            {
                est_taglio[pos]=mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[lato][0]][1];
                pos++;
            }
        }
        //----------------------ora aggiorno le celle 1D e infine anche le celle 2D----------------------
        //celle 1D
        unsigned int ID_1D=mesh.NumberCell1D;
        mesh.Cell1DId.push_back(ID_1D);
        mesh.NumberCell1D++;
        mesh.Cell1DMarkers.push_back(0); //il marker sicuramente è nullo dato che il segmento taglia la frattura (il caso degenere è stato già trattato prima)
        mesh.Cell1DVertices.push_back(est_taglio);

        //celle 2D
        mesh.Cell2DId.push_back(mesh.NumberCell2D);
        mesh.NumberCell2D++;

        //----------------------creazione dei 2 nuovi poligoni a partire da quello vecchio----------------------
        //----------------------POLIGONO 1----------------------

        // (1) inizializzo il poligono con tutto quello che so
        std::vector<unsigned int> nuovo_poligono;
        unsigned int larghezza=mesh.Cell2DEdges[poligoni[i]].size();
        nuovo_poligono.reserve(larghezza);
        nuovo_poligono.push_back(ID_1D); //come primo lato metto il nuovo lato che taglia il vecchio poligono
        //dei 2 lati del vecchio poligono tagliati, scelgo di andare sul primo
        const unsigned int& vecchio_lato_1=lati_coinvolti[i][0];
        const unsigned int& vecchio_lato_2=lati_coinvolti[i][1];

        // (2) capisco come partire per visitare il nuovo poligono in senso antiorario come si faceva nel vecchio poligono
        //trovo dove in quale posizione di mesh.Cell2DEdges[poligoni[i]] si trova vecchio_lato_1
        auto it=std::find(mesh.Cell2DEdges[poligoni[i]].begin(),mesh.Cell2DEdges[poligoni[i]].end(),vecchio_lato_1);
        unsigned int posto=std::distance(mesh.Cell2DEdges[poligoni[i]].begin(),it);
        posto++;
        posto=posto%larghezza;
        unsigned int lato_per_test=mesh.Cell2DEdges[poligoni[i]][posto];

        //ste info servono per velocizzare il processo
        int scelta_partenza_lato_nuovo; //posizione nel vettore dei nuovi lati che sostituiscono il vecchio lato
        int scelta_arrivo_lato_nuovo; //posizione nel vettore dei nuovi lati che sostituiscono il vecchio lato
        bool fine;

        //test effettivo per capire come partire: ho tanti casi da tenere in considerazione
        //CREDO SI POSSA TOGLIERE QUALCOSA
        //LPT=D?
        if(lato_per_test==vecchio_lato_2)
        { //S
            //A-D?
            if(!(std::find(mesh.Cell1DVertices[lato_per_test].begin(),mesh.Cell1DVertices[lato_per_test].end(),mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_1][0]][0])==mesh.Cell1DVertices[lato_per_test].end()) || !(std::find(mesh.Cell1DVertices[lato_per_test].begin(),mesh.Cell1DVertices[lato_per_test].end(),mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_1][0]][1])==mesh.Cell1DVertices[lato_per_test].end()))
            { //S
                //+A
                nuovo_poligono.push_back(mappa_vecchi_lati_nuovi_lati[vecchio_lato_1][0]);
                //+D
                nuovo_poligono.push_back(lato_per_test);
                //ho finito
                scelta_partenza_lato_nuovo=-1;
                scelta_arrivo_lato_nuovo=-1;
                fine=true;
            }
            else
            { //N
                //B-D?
                if(!(std::find(mesh.Cell1DVertices[lato_per_test].begin(),mesh.Cell1DVertices[lato_per_test].end(),mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_1][1]][0])==mesh.Cell1DVertices[lato_per_test].end()) || !(std::find(mesh.Cell1DVertices[lato_per_test].begin(),mesh.Cell1DVertices[lato_per_test].end(),mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_1][1]][1])==mesh.Cell1DVertices[lato_per_test].end()))
                { //S
                    //+B
                    nuovo_poligono.push_back(mappa_vecchi_lati_nuovi_lati[vecchio_lato_1][1]);
                    //+D
                    nuovo_poligono.push_back(lato_per_test);
                    //ho finito
                    scelta_partenza_lato_nuovo=-1;
                    scelta_arrivo_lato_nuovo=0;
                    fine=true;
                }
                else
                { //N
                    //A-C?
                    if(!(std::find(mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_2][1]].begin(),mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_2][1]].end(),mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_1][0]][0])==mesh.Cell1DVertices[lato_per_test].end()) || !(std::find(mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_2][1]].begin(),mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_2][1]].end(),mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_1][0]][1])==mesh.Cell1DVertices[lato_per_test].end()))
                    { //S
                        //+A
                        nuovo_poligono.push_back(mappa_vecchi_lati_nuovi_lati[vecchio_lato_1][0]);
                        //+C
                        nuovo_poligono.push_back(mappa_vecchi_lati_nuovi_lati[vecchio_lato_2][1]);
                        //ho finito
                        scelta_partenza_lato_nuovo=1;
                        scelta_arrivo_lato_nuovo=0;
                        fine=true;
                    }
                    else
                    { //N
                        //+B
                        nuovo_poligono.push_back(mappa_vecchi_lati_nuovi_lati[vecchio_lato_1][1]);
                        //+C
                        nuovo_poligono.push_back(mappa_vecchi_lati_nuovi_lati[vecchio_lato_2][1]);
                        //ho finito
                        scelta_partenza_lato_nuovo=0;
                        scelta_arrivo_lato_nuovo=0;
                        fine=true;
                    }
                }
            }
        }
        else
        { //N
            //A-LPT?
            if(!(std::find(mesh.Cell1DVertices[lato_per_test].begin(),mesh.Cell1DVertices[lato_per_test].end(),mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_1][0]][0])==mesh.Cell1DVertices[lato_per_test].end()) || !(std::find(mesh.Cell1DVertices[lato_per_test].begin(),mesh.Cell1DVertices[lato_per_test].end(),mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_1][0]][1])==mesh.Cell1DVertices[lato_per_test].end()))
            { //S
                //T-LPT?
                if(!(std::find(mesh.Cell1DVertices[lato_per_test].begin(),mesh.Cell1DVertices[lato_per_test].end(),est_taglio[0])==mesh.Cell1DVertices[lato_per_test].end()) || !(std::find(mesh.Cell1DVertices[lato_per_test].begin(),mesh.Cell1DVertices[lato_per_test].end(),est_taglio[1])==mesh.Cell1DVertices[lato_per_test].end()))
                { //S
                    //+LPT
                    nuovo_poligono.push_back(lato_per_test);
                    posto++;
                    posto=posto%larghezza;
                    scelta_arrivo_lato_nuovo=0;
                }
                else
                { //N
                    //+A
                    nuovo_poligono.push_back(mappa_vecchi_lati_nuovi_lati[vecchio_lato_1][0]);
                    //+LPT
                    nuovo_poligono.push_back(lato_per_test);
                    posto++;
                    posto=posto%larghezza;
                    scelta_arrivo_lato_nuovo=-1;
                }
            }
            else
            { //N
                //+B
                nuovo_poligono.push_back(mappa_vecchi_lati_nuovi_lati[vecchio_lato_1][1]);
                //+LPT
                nuovo_poligono.push_back(lato_per_test);
                posto++;
                posto=posto%larghezza;
                scelta_arrivo_lato_nuovo=0;
            }
        }

        if(!fine) //se non ho finito
        {
            // (3) inserisco i nuovi lati finchè non trovo l'altro lato vecchio coinvolto dal taglio
            while(std::find(lati_coinvolti[i].begin(),lati_coinvolti[i].end(),mesh.Cell2DEdges[poligoni[i]][posto])==lati_coinvolti[i].end())
            {
                nuovo_poligono.push_back(mesh.Cell2DEdges[poligoni[i]][posto]);
                posto++;
                posto=posto%larghezza;
            }
            lato_per_test=nuovo_poligono.back();

            // (4) capisco come finire
            //LPT-D?
            if(!(std::find(mesh.Cell1DVertices[lato_per_test].begin(),mesh.Cell1DVertices[lato_per_test].end(),mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_2][0]][0])==mesh.Cell1DVertices[lato_per_test].end()) || !(std::find(mesh.Cell1DVertices[lato_per_test].begin(),mesh.Cell1DVertices[lato_per_test].end(),mesh.Cell1DVertices[mappa_vecchi_lati_nuovi_lati[vecchio_lato_2][0]][1])==mesh.Cell1DVertices[lato_per_test].end()))
            { //S
                //T-/-LPT?
                if(std::find(mesh.Cell1DVertices[lato_per_test].begin(),mesh.Cell1DVertices[lato_per_test].end(),est_taglio[0])==mesh.Cell1DVertices[lato_per_test].end() && std::find(mesh.Cell1DVertices[lato_per_test].begin(),mesh.Cell1DVertices[lato_per_test].end(),est_taglio[1])==mesh.Cell1DVertices[lato_per_test].end())
                { //S
                    //+D
                    nuovo_poligono.push_back(mappa_vecchi_lati_nuovi_lati[vecchio_lato_2][0]);
                    scelta_partenza_lato_nuovo=-1;
                }
                else
                { //N
                    scelta_partenza_lato_nuovo=0;
                }
            }
            else
            { //N
                //+C
                nuovo_poligono.push_back(mappa_vecchi_lati_nuovi_lati[vecchio_lato_2][1]);
                scelta_partenza_lato_nuovo=0;
            }
        }

        // (5) salvo il nuovo poligono fatto da lati nella mesh
        mesh.Cell2DEdges.push_back(nuovo_poligono);

        // (6) ho il poligono descritto come un vettore di lati, ora costruisco il poligono come vettore di vertici a partire da quello come vettore di lati
        std::vector<unsigned int> vertici_nuovo_poligono;
        vertici_nuovo_poligono.reserve(larghezza);
        for(unsigned int i=0;i<nuovo_poligono.size();i++)
        {
            if(!(vertici_nuovo_poligono.empty()))
            {
                if(mesh.Cell1DVertices[nuovo_poligono[i-1]][1]==mesh.Cell1DVertices[nuovo_poligono[i]][0])
                {
                    vertici_nuovo_poligono.push_back(mesh.Cell1DVertices[nuovo_poligono[i]][0]);
                }
                else
                {
                    if(mesh.Cell1DVertices[nuovo_poligono[i-1]][0]==mesh.Cell1DVertices[nuovo_poligono[i]][0])
                    {
                        vertici_nuovo_poligono.push_back(mesh.Cell1DVertices[nuovo_poligono[i]][0]);
                    }
                    else
                    {
                        vertici_nuovo_poligono.push_back(mesh.Cell1DVertices[nuovo_poligono[i]][1]);
                    }
                }
            }
            else
            {
                //FORSE TUTTO STO COSO NON SERVE
                if(mesh.Cell1DVertices[nuovo_poligono[(i+1)%larghezza]][0]==mesh.Cell1DVertices[nuovo_poligono[i]][1] || mesh.Cell1DVertices[nuovo_poligono[(i+1)%larghezza]][1]==mesh.Cell1DVertices[nuovo_poligono[i]][1])
                {
                    vertici_nuovo_poligono.push_back(mesh.Cell1DVertices[nuovo_poligono[i]][0]);
                }
                else
                {
                    vertici_nuovo_poligono.push_back(mesh.Cell1DVertices[nuovo_poligono[i]][1]);
                }
            }
        }

        // (7) salvo il poligono fatto da vertici della mesh
        mesh.Cell2DVertices.push_back(vertici_nuovo_poligono);

        //----------------------POLIGONO 2----------------------
        //per risparmiare memoria, uso la struttura dati del poligono 1
        nuovo_poligono.clear();
        nuovo_poligono.push_back(ID_1D);
        if(scelta_partenza_lato_nuovo!=-1)
        { //S
            nuovo_poligono.push_back(mappa_vecchi_lati_nuovi_lati[vecchio_lato_2][scelta_partenza_lato_nuovo]);
        }
        else
        { //N
            //D=D,C?
            if(mappa_vecchi_lati_nuovi_lati[vecchio_lato_2].size()==2)
            { //S
                nuovo_poligono.push_back(mappa_vecchi_lati_nuovi_lati[vecchio_lato_2][1]);
            }
        }
        it=std::find(mesh.Cell2DEdges[poligoni[i]].begin(),mesh.Cell2DEdges[poligoni[i]].end(),vecchio_lato_2);
        posto=std::distance(mesh.Cell2DEdges[poligoni[i]].begin(),it);
        posto++;
        posto=posto%larghezza;
        while(std::find(lati_coinvolti[i].begin(),lati_coinvolti[i].end(),mesh.Cell2DEdges[poligoni[i]][posto])==lati_coinvolti[i].end())
        {
            nuovo_poligono.push_back(mesh.Cell2DEdges[poligoni[i]][posto]);
            posto++;
            posto=posto%larghezza;
        }
        if(scelta_arrivo_lato_nuovo!=-1)
        { //S
            nuovo_poligono.push_back(mappa_vecchi_lati_nuovi_lati[vecchio_lato_1][scelta_arrivo_lato_nuovo]);
        }
        else
        { //N
            //A=A,B?
            if(mappa_vecchi_lati_nuovi_lati[vecchio_lato_1].size()==2)
            { //S
                nuovo_poligono.push_back(mappa_vecchi_lati_nuovi_lati[vecchio_lato_1][1]);
            }
        }
        mesh.Cell2DEdges[poligoni[i]]=nuovo_poligono;

        vertici_nuovo_poligono.clear();
        for(unsigned int i=0;i<nuovo_poligono.size();i++)
        {
            if(!(vertici_nuovo_poligono.empty()))
            {
                if(mesh.Cell1DVertices[nuovo_poligono[i-1]][1]==mesh.Cell1DVertices[nuovo_poligono[i]][0])
                {
                    vertici_nuovo_poligono.push_back(mesh.Cell1DVertices[nuovo_poligono[i]][0]);
                }
                else
                {
                    if(mesh.Cell1DVertices[nuovo_poligono[i-1]][0]==mesh.Cell1DVertices[nuovo_poligono[i]][0])
                    {
                        vertici_nuovo_poligono.push_back(mesh.Cell1DVertices[nuovo_poligono[i]][0]);
                    }
                    else
                    {
                        vertici_nuovo_poligono.push_back(mesh.Cell1DVertices[nuovo_poligono[i]][1]);
                    }
                }
            }
            else
            {
                //FORSE TUTTO STO COSO NON SERVE
                if(mesh.Cell1DVertices[nuovo_poligono[(i+1)%larghezza]][0]==mesh.Cell1DVertices[nuovo_poligono[i]][1] || mesh.Cell1DVertices[nuovo_poligono[(i+1)%larghezza]][1]==mesh.Cell1DVertices[nuovo_poligono[i]][1])
                {
                    vertici_nuovo_poligono.push_back(mesh.Cell1DVertices[nuovo_poligono[i]][0]);
                }
                else
                {
                    vertici_nuovo_poligono.push_back(mesh.Cell1DVertices[nuovo_poligono[i]][1]);
                }
            }
        }
        mesh.Cell2DVertices[poligoni[i]]=vertici_nuovo_poligono;
    }
}




//funzione che trova i poligoni che intersecano il segmento
//OUTPUT: tupla di 2 elementi: 1) vettore degli indici dei poligoni intersecati dal segmento 2)vettori di vettori dei lati che intersecano l'asse che contiene il segmento
std::tuple<std::vector<unsigned int>,std::vector<std::array<unsigned int,2>>> contatto_poligoni_segmento(const DFN& disc_frac_net,const unsigned int& id,PolygonalMesh& mesh,const Vector3d& est1,const Vector3d& est2, const double& tol)
{

    //riga 1 e riga 3 della matrice di rotazione per spostare il problema sul piano yz (la riga 2 non serve)
    Vector3d vec1_rot=(est2-est1);
    double norma=vec1_rot.norm();
    vec1_rot.normalize();
    Vector3d vec2_rot;
    vec2_rot=(vec1_rot.cross(disc_frac_net.versori[id]));

    std::vector<unsigned int> temp_pol; //output 1
    std::vector<std::array<unsigned int,2>> temp_lat; //output 2
    temp_pol.reserve(mesh.NumberCell2D);
    temp_lat.reserve(mesh.NumberCell2D);

    //selezione dei poligoni che intersecano il segmento
    for(unsigned int n2d=0; n2d < mesh.NumberCell2D; n2d++) //ciclo su tutti i poligoni della mesh
    {
        //vettore dei lati che intersecano l'asse del segmento
        std::array<unsigned int,2> lati_coinvolti={std::numeric_limits<unsigned int>::max(),std::numeric_limits<unsigned int>::max()}; //è improbabile che ci sia veramente un lato coinvolto che abbia questo indice
        unsigned int pos=0;

        //estremi del primo lato del poligono
        double x_1=vec1_rot.dot(mesh.Cell0DCoordinates[mesh.Cell2DVertices[n2d][0]]-est1);        
        double z_1=vec2_rot.dot(mesh.Cell0DCoordinates[mesh.Cell2DVertices[n2d][0]]-est1);

        double x_2=vec1_rot.dot(mesh.Cell0DCoordinates[mesh.Cell2DVertices[n2d][1]]-est1);
        double z_2=vec2_rot.dot(mesh.Cell0DCoordinates[mesh.Cell2DVertices[n2d][1]]-est1);

        //estremi della sezione
        double min_sez=INFINITY;
        double max_sez=-INFINITY;

        double intersezione_1=NAN; //eventuale punto di intersezione del lato con la retta
        if(!(z_2*z_1>tol)) //se il lato interseca la retta che contiene il segmento
        {
            lati_coinvolti[pos]=mesh.Cell2DEdges[n2d][0];
            pos++;
            intersezione_1=x_2-z_2*(x_1-x_2)/(z_1-z_2); //punto d'intersezione del lato con la retta
            min_sez=intersezione_1;
            max_sez=intersezione_1;
        }
        //riinizializzo x_1 e z_1 per passare al lato successivo
        x_1=x_2;
        z_1=z_2;
        for(unsigned int id_vertice=1;id_vertice<mesh.Cell2DVertices[n2d].size();id_vertice++) //ciclo su tutti i lati restanti del poligono
        {
            //stessa roba di prima
            x_2=vec1_rot.dot(mesh.Cell0DCoordinates[mesh.Cell2DVertices[n2d][(id_vertice+1)%mesh.Cell2DVertices[n2d].size()]]-est1);
            z_2=vec2_rot.dot(mesh.Cell0DCoordinates[mesh.Cell2DVertices[n2d][(id_vertice+1)%mesh.Cell2DVertices[n2d].size()]]-est1);
            double intersezione_2=NAN;
            if(!(z_2*z_1>tol))
            {
                intersezione_2=x_2-z_2*(x_1-x_2)/(z_1-z_2);
                if(!std::isnan(intersezione_1)) //se intersezione_1 non è NAN faccio il controllo dopo
                {
                    if(std::abs(intersezione_2-intersezione_1)/std::max(std::max(intersezione_2,intersezione_1),1.)>tol && pos!=2) //sto controllo serve quando ho l'intersezione della retta con un lato che avviene proprio in un estremo del lato e va così a coinvolgere anche il lato successivo
                    {
                        lati_coinvolti[pos]=mesh.Cell2DEdges[n2d][id_vertice];
                        pos++;
                        min_sez=std::min(min_sez,intersezione_2);
                        max_sez=std::max(max_sez,intersezione_2);
                    }
                }
                else //se invece intersezione_1 è NAN non serve fare il controllo
                {
                    lati_coinvolti[pos]=mesh.Cell2DEdges[n2d][id_vertice];
                    pos++;
                    min_sez=std::min(min_sez,intersezione_2);
                    max_sez=std::max(max_sez,intersezione_2);
                }
            }
            x_1=x_2;
            z_1=z_2;
            intersezione_1=intersezione_2;
        }

        //se il poligono interseca il segmento
        if(!(max_sez<tol || min_sez/std::max(min_sez,norma)>norma/std::max(min_sez,norma)-tol) && std::abs(max_sez-min_sez)/std::max(max_sez,min_sez)>tol)
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
        //salvo le informazioni che so già sulla mesh
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

            pol_vertici_lati.push_back(i);
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
            std::vector<unsigned int> poligoni_coinvolti=std::get<std::vector<unsigned int>>(poligoni_e_lati_coinvolti);
            std::vector<std::array<unsigned int,2>> lati_coinvolti=std::get<std::vector<std::array<unsigned int,2>>>(poligoni_e_lati_coinvolti);
            aggiorna_mesh(mesh,poligoni_coinvolti,lati_coinvolti,est1,est2,tol);
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
