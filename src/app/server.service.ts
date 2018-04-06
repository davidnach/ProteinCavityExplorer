import {Injectable} from '@angular/core';
import {HttpClient} from '@angular/common/http';

@Injectable()
export class ServerService{
  constructor(private http: HttpClient){}


    createExperiment(pid : {pid :string}){
      return this.http.post('http://localhost:3000/api/createExperiment', pid);
    }

    checkIfPidExists(pid: string){
      return this.http.get('https://files.rcsb.org/view/1ake.pdb')
    }

}
