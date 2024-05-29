#include"Utils.hpp"
#include <cassert>

using namespace Eigen;

namespace LibraryDFN{




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




//funzione che scarta le fratture che sicuramente non si intersecano, l'output è un vettore di indici degli indici di DFN.idfratture che non sono stati scartati, quindi su cui bisogna controllare manualmente se si intersechino
std::vector<std::array<unsigned int,2>> scarta_fratture(DFN& disc_frac_net)
{
    std::vector<std::array<unsigned int,2>> indici_validi; //vettore finale di output con gli indici che non hanno passato il test
    //ATTENZIONE ALL'OVERFLOW NEL PROSSIMO PASSAGGIO, PROBABILMENTE E' DA CAMBIARE
    indici_validi.reserve(disc_frac_net.numFratture*(disc_frac_net.numFratture-1)/2); //al massimo ho n su 2 coppie
    std::vector<std::tuple<unsigned int,Vector3d,double>> info_baricentro; //vettore con tutte le info da confrontare per effettuare il test
    info_baricentro.reserve(disc_frac_net.numFratture);

    for(unsigned int i=0;i<disc_frac_net.numFratture;i++) // in questo ciclo salvo le info per il confronto
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

    //test effettivo per scartare le fratture che sicuramente non si intersecano
    for(unsigned int i=0;i<disc_frac_net.numFratture-1;i++) //con questi 2 cicli for prendo tutte le combinazioni di 2 elementi delle fratture per il confronto, senza ripetizioni
    {
        for(unsigned int j=i+1;j<disc_frac_net.numFratture;j++)
        {
            std::array<unsigned int,2> coppia={i,j};
            if((std::get<Vector3d>(info_baricentro[i])-std::get<Vector3d>(info_baricentro[j])).squaredNorm()<std::get<double>(info_baricentro[i])+std::get<double>(info_baricentro[j]))
            {
                indici_validi.push_back(coppia);
            }
        }
    }
    return indici_validi;
}




//funzione che genera la matrice di rotazione che trasforma il piano in cui la frattura è contenuta nel piano xy, in modo che il primo lato sia parallelo all'asse x (è una scelta, se ne poteva fare anche un'altra)
Eigen::Matrix<double,3,3> allinea_xy(const std::vector<Eigen::Vector3d>& vertici)
{
    Eigen::Matrix<double,3,3> rotazione;
    Eigen::Vector3d n=(vertici[0]-vertici[1]).cross(vertici[2]-vertici[1]); //forse è un po' grossolana come stima del versore normale, tuttavia serve solamente per rendere più facile verificare una condizione e ridurre il costo computazionale
    n=1./n.norm()*n;
    rotazione.row(2)=n;
    Eigen::Vector3d v=vertici[0]-vertici[1];
    v=1/v.norm()*v;
    rotazione.row(0)=v;
    rotazione.row(1)=n.cross(v);
    return rotazione;
}




//funzione effettiva che memorizza le tracce nella struttura DFN
void memorizza_tracce(DFN& disc_frac_net,double tol)
{
    //controllo che la tolleranza scelta dall'utente non sia minore della precisione di macchina
    tol=std::max(tol,__DBL_EPSILON__);

    //--------- passo 1 ---------
    //scrematura grossolana delle fratture che sicuramente non si intersecano
    std::vector<std::array<unsigned int,2>> fratture_suscettibili=scarta_fratture(disc_frac_net);

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
        std::array<Eigen::Vector3d,2> est_tra={Eigen::Vector3d(NAN,NAN,NAN),Eigen::Vector3d(NAN,NAN,NAN)};
        unsigned int pos=0; //è l'indice di riempimento di est_tra, ossia indica la posizione libera di est_tra (se pos==2 vuol dire che è pieno)

        //--------- passo 4 ---------
        //ruoto i poligoni sul piano xy per rendere più facili i confronti dopo
        Eigen::Matrix<double,3,3> mat_rot=allinea_xy(disc_frac_net.vertici[coppia[0]]);
        //std::vector<Eigen::Vector3d> vertici_ruotati_1={Eigen::Vector3d(0.,0.,0.)}; //prima frattura
        std::vector<Eigen::Vector3d> vertici_ruotati_2={}; //seconda frattura
        vertici_ruotati_2.resize(disc_frac_net.numVertici[coppia[1]]);

        // for(unsigned int i=1;i<disc_frac_net.numVertici[coppia[0]];i++)
        // {
        //     vertici_ruotati_1[i]=mat_rot*(disc_frac_net.vertici[coppia[0]][i]-disc_frac_net.vertici[coppia[0]][0]);
        // }
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
                    //std::cout << "nuovo triangolo analizzato"<<std::endl;
                    Eigen::Matrix<double,3,3> mat_coeff;
                    Eigen::Vector3d b=disc_frac_net.vertici[coppia[1]][i]-disc_frac_net.vertici[coppia[0]][triangolo[0]]; //vettore dei termini noti
                    mat_coeff.col(0)=disc_frac_net.vertici[coppia[0]][triangolo[1]]-disc_frac_net.vertici[coppia[0]][triangolo[0]];
                    mat_coeff.col(1)=disc_frac_net.vertici[coppia[0]][triangolo[2]]-disc_frac_net.vertici[coppia[0]][triangolo[0]];
                    mat_coeff.col(2)=disc_frac_net.vertici[coppia[1]][i]-disc_frac_net.vertici[coppia[1]][(i+1)%disc_frac_net.numVertici[coppia[1]]];
                    Eigen::Vector3d x=mat_coeff.lu().solve(b);
                    //std::cout << x<<std::endl;

                    //--------- passo 8 ---------
                    //controllo se la soluzione trovata appartiene al triangolo
                    if(x(0)>=-tol && x(1)>=-tol && x(0)+x(1)<=1.+tol && x(2)<=1.+tol && x(2)>=-tol) //l'ultimo controllo non dovrebbe servire ma per ora lo tengo
                    {
                        Eigen::Vector3d y=0.5*(disc_frac_net.vertici[coppia[1]][i]+x(2)*(-mat_coeff.col(2))+disc_frac_net.vertici[coppia[0]][triangolo[0]]+x(0)*mat_coeff.col(0)+x(1)*mat_coeff.col(1));
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
            mat_rot=allinea_xy(disc_frac_net.vertici[coppia[1]]);
            std::vector<Eigen::Vector3d> vertici_ruotati_1;
            vertici_ruotati_1.resize(disc_frac_net.numVertici[coppia[0]]);
            // for(unsigned int i=1;i<disc_frac_net.numVertici[coppia[1]];i++)
            // {
            //     vertici_ruotati_2[i]=mat_rot*(disc_frac_net.vertici[coppia[1]][i]-disc_frac_net.vertici[coppia[1]][0]);
            // }
            for(unsigned int i=0;i<disc_frac_net.numVertici[coppia[0]];i++)
            {
                vertici_ruotati_1[i]=mat_rot*(disc_frac_net.vertici[coppia[0]][i]-disc_frac_net.vertici[coppia[1]][0]);
            }

            triangolazione=triangola_frattura(disc_frac_net,coppia[1]);

            for(unsigned int i=0;i<disc_frac_net.numVertici[coppia[0]];i++)
            {
                //std::cout << i <<std::endl;
                if(!(vertici_ruotati_1[i](2)*vertici_ruotati_1[(i+1)%disc_frac_net.numVertici[coppia[0]]](2)>tol))
                {
                    for(auto& triangolo:triangolazione)
                    {
                        //std::cout << "nuovo triangolo analizzato"<<std::endl;
                        Eigen::Matrix<double,3,3> mat_coeff;
                        Eigen::Vector3d b=disc_frac_net.vertici[coppia[0]][i]-disc_frac_net.vertici[coppia[1]][triangolo[0]];
                        mat_coeff.col(0)=disc_frac_net.vertici[coppia[1]][triangolo[1]]-disc_frac_net.vertici[coppia[1]][triangolo[0]];
                        mat_coeff.col(1)=disc_frac_net.vertici[coppia[1]][triangolo[2]]-disc_frac_net.vertici[coppia[1]][triangolo[0]];
                        mat_coeff.col(2)=disc_frac_net.vertici[coppia[0]][i]-disc_frac_net.vertici[coppia[0]][(i+1)%disc_frac_net.numVertici[coppia[0]]];
                        Eigen::Vector3d x=mat_coeff.lu().solve(b);
                        //std::cout << x<<std::endl;
                        if(x(0)>=-tol && x(1)>=-tol && x(0)+x(1)<=1.+tol && x(2)<=1.+tol && x(2)>=-tol)
                        {
                            Eigen::Vector3d y=0.5*(disc_frac_net.vertici[coppia[0]][i]+x(2)*(-mat_coeff.col(2))+disc_frac_net.vertici[coppia[1]][triangolo[0]]+x(0)*mat_coeff.col(0)+x(1)*mat_coeff.col(1));
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
            disc_frac_net.tracce.push_back(coppia);
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
        //std::cout<<"++++++++++++++++++++++++"<<std::endl;
    }
    disc_frac_net.numTracce=idt;
    std::cout << "Alla fine sono state trovate "<<idt<<" tracce"<<std::endl;
}
//=====================================================================================================




//funzione che trova l'intersezione tra 2 segmenti
std::tuple<Vector3d,bool> interseca_segmenti(const Vector3d& A1,const Vector3d& A2,const Vector3d& B1,const Vector3d& B2,const double& tol)
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
    bool interno=(alfa<1+tol && alfa>-tol && beta<1+tol && beta>-tol);
    return std::make_tuple(x,interno);
}
}
