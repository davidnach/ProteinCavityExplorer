import { Component } from '@angular/core';
import {Router} from '@angular/router';
import {ServerService} from './server.service';
import {UserInputService} from './userinput.service';
import {GlobalData} from './globaldata.service';

@Component({
  selector: 'app-root',
  templateUrl: './app.component.html',
  styleUrls: ['./app.component.css']
})
export class AppComponent {
  constructor(private serverService: ServerService,
  private router: Router,private userInputService: UserInputService,
private globaldata : GlobalData ){

    this.userInputService.pidEntered.subscribe(
       (pid : string) => {
         this.submitPidToServer(pid);
          });

   this.userInputService.experimentIdEntered.subscribe(
	(expId : string) => {
	  this.submitExperimentId(expId);
	});
   }

  onSubmit(){
    this.serverService.createExperiment({pid: 'lake'})
    .subscribe((response) => {
      console.log(response);
      console.log(response['pid'])},
    (error) => {console.log(error)});
  }

  submitPidToServer(pid : string){
    this.serverService.createExperiment({pid: 'lake'})
    .subscribe((response) => {
      console.log(response['expId']);
      this.globaldata.setExperimentId(response['expId']);
	this.router.navigate(['experiment']);},
    (error) => {console.log(error)});
  }

  submitExperimentId(experimentId : string){
    this.serverService.retrieveExperiment({expId : '11111'})
    .subscribe((response) => {
     console.log(response);
     this.globaldata.setExperimentId(response['expId']);},
    (error) => {console.log(error)});
    this.router.navigate(['experiment']);
  }



}
