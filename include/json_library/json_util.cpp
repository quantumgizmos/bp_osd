//
// Created by joschka on 30/04/2020.
//

#include "json_util.h"

json json_read_safe(json j, string key, string default_to){

    try {
        auto data_item = j[key];

        if(data_item==nullptr) {
            if(default_to=="0") throw 22;
            cout<<"Key: "<<key<< " not found. Defaulting to: "<<default_to<<" method"<<endl;
            return default_to;
        }

        return data_item;
    }
    catch (int e)
    {
        // output exception information
        cout<<"Error. Key '"<<key<<"' does not exist!"<<endl;
        exit(23);
    }

}

