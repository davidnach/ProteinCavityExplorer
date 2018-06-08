/*all communications with the server go here*/


import {Injectable} from '@angular/core';
import {HttpClient} from '@angular/common/http';

@Injectable()

export class ServerService{
  constructor(private http: HttpClient){}


    createExperiment(pid : {pid :string}){
      return this.http.post('/api/createExperiment', pid);
    }

   
    retrieveExperiment(experimentId : {expId : string}){
    	return this.http.post('/api/retrieveExperiment',experimentId);
    }


     getPidsWithSurfReports(){
	
	return this.http.post('/api/getPidsWithSurfReports',{ignore : "ignore"});
    }

    
}
