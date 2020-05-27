//
// Created by joschka on 30/04/2020.
//

#include "json_util.h"

json json_read_safe(json j, string key){

    try {
        auto data_item = j[key];

        if(data_item==nullptr) {
            throw 22;
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

