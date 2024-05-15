#include"Utils.hpp"

using namespace LibraryDFN;
using namespace Eigen;

//-----------------------------------------------------------------------------------------
//qua si triangola la singola frattura (sticazzi se vengono triangoli degeneri, anche se non dovrebbero venire)
std::vector<std::array<unsigned int,3>> triangola_frattura(DFN& disc_frac_net,unsigned int idr) //idr=id relativo, nel senso che è l'indice del vettore numvertices
{
    std::vector<std::array<unsigned int,3>> triangolazione;
    triangolazione.reserve(disc_frac_net.numVertici[idr]-2);
    for(unsigned int i=1;i<disc_frac_net.numVertici[idr]-1;i++)
    {
        std::array<unsigned int,3> triangolo={0,i,i+1};
        triangolazione.push_back(triangolo);
    }
    return triangolazione
}
//-----------------------------------------------------------------------------------------
//QUESTO CODICE QUA SOTTO VA OTTIMIZZATO
//qua si scartano gli id delle fratture che sicuramente non si intersecano: l'output è un vettore di indici degli indici di DFN.idfratture che non sono stati scartati, quindi su cui bisogna controllare manualmente se si intersechino
std::vector<std::array<unsigned int,2>> scarta_fratture(DFN& disc_frac_net)
{
    std::vector<std::array<unsigned int,2>> indici_validi; //vettore finale di output con gli indici che non hanno passato il test
    std::list<std::array<unsigned int,2>> temporanea;
    std::vector<std::tuple<unsigned int,Vector3d,double>> info_baricentro; //vettore con tutte le info da confrontare per effettuare il test
    info_baricentro.reserve(disc_frac_net.numFratture);

    for(unsigned int i=0;i<disc_frac_net.numFratture;i++) //salvo le info per il confronto
    {
        std::cout << "poligono: "<< i<<std::endl;
        Vector3d baricentro(0,0,0);
        for(Vector3d& coord:disc_frac_net.vertici[i])
        {
            baricentro+=coord;
        }
        baricentro=baricentro/disc_frac_net.numVertici[i];
        //std::cout << " baricentro: "<<std::endl<< baricentro<<std::endl;
        double raggio=0.;
        for(Vector3d& coord:disc_frac_net.vertici[i])
        {
            double candidato_raggio=(coord-baricentro).norm();
            if(candidato_raggio>raggio)
                raggio=candidato_raggio;
        }
        //std::cout << " raggio: "<< raggio << std::endl;
        std::tuple<unsigned int,Vector3d,double> tupla(i,baricentro,raggio);
        info_baricentro.push_back(tupla);
    }

    //test effettivo per scartare le fratture che sicuramente non si intersecano
    for(unsigned int i=0;i<disc_frac_net.numFratture-1;i++) //con questi 2 cicli for prendo tutte le combinazioni di 2 elementi delle fratture per il confronto, senza ripetizioni
    {
        for(unsigned int j=i+1;j<disc_frac_net.numFratture;j++)
        {
            //std::cout<<"("<<i<<","<<j<<")"<<std::endl;
            std::array<unsigned int,2> coppia={i,j};
            if((std::get<Vector3d>(info_baricentro[i])-std::get<Vector3d>(info_baricentro[j])).norm()<std::get<double>(info_baricentro[i])+std::get<double>(info_baricentro[j])) //forse ci va l'= anche
            {
                temporanea.push_back(coppia);
                //std::cout<<"("<<coppia[0]<<","<<coppia[1]<<")"<<std::endl;
            }
        }
    }
    //CAMBIARE STA COSA ASSOLUTAMENTE
    indici_validi.reserve(temporanea.size());
    for (auto& coppia : temporanea)
    {
        indici_validi.push_back(coppia);
    }
    return indici_validi;
}
//-----------------------------------------------------------------------------------------
//funzione effettiva che memorizza le tracce
/*void memorizza_tracce(DFN& disc_frac_net)
{
    //scrematura delle fratture che sicuramente non si intersecano
    std::vector<unsigned int> fratture_suscettibili;
    fratture_suscettibili=scarta_fratture(disc_frac_net);

    //ricerca degli estremi della traccia
    for(auto& coppia:fratture_suscettibili)
    {
        std::vector<std::array<unsigned int,3>> triangolazione=triangola_frattura(disc_frac_net,coppia[0]);
        for(auto& triangolo:triangolazione)
        {

        }
    }

}*/
