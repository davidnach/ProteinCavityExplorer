/*user-input component validation and commnication w/ other components*/

import { Injectable,EventEmitter} from '@angular/core';
import {HttpClient} from '@angular/common/http';

@Injectable()

export class UserInputService{

  pid:string;

  pidEntered = new EventEmitter<string>();
  experimentIdEntered = new EventEmitter<string>();

  constructor(private http:HttpClient){ }

  checkIfPidExists(pid: string){
  //getPidStatus(pid: string){
    //return this.http.get('https://www.rcsb.org/pdb/rest/idStatus?structureId='+pid+'');
      return this.http.get('https://files.rcsb.org/view/' + pid + '.pdb');
      //return this.http.get('https://files.rcsb.org/view/1ake.pdb');
  }

  setPid(pid : string){
    this.pid = pid;
    console.log(this.pid);
  }

  getPid(){
    return this.pid;
  }


}
