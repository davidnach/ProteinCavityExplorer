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

  submitPidToServer(pid : string){
   
    this.serverService.createExperiment({pid: pid})
    .subscribe((response) => {
      console.log(response);
      console.log(response['expId']);
      console.log(response['pockets']);
      this.globaldata.setExperimentId(response['expId']);
      this.globaldata.setPockets(response['pockets']);	
	this.router.navigate(['experiment']);},
    (error) => {console.log(error)});
  
   this.router.navigate(['experiment']);

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
