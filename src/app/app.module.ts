import { BrowserModule } from '@angular/platform-browser';
import { NgModule } from '@angular/core';
import {HttpClientModule} from '@angular/common/http';
import { FormsModule } from '@angular/forms';
import {ServerService} from './server.service';
import {UserInputService} from './userinput.service';
import {Routes, RouterModule} from "@angular/router"


import { AppComponent } from './app.component';
import { UserInputComponent } from './user-input/user-input.component';
import { ExperimentResultsComponent} from './experiment-results/experiment-results.component';

const appRoutes: Routes = [

  {path: '', redirectTo: 'home', pathMatch: 'full'},
  {path: 'home',component: UserInputComponent},
  {path: 'experiment', component: ExperimentResultsComponent}
];

@NgModule({
  declarations: [
    AppComponent,
    UserInputComponent,
    ExperimentResultsComponent
  ],
  imports: [
    BrowserModule,
    HttpClientModule,
    FormsModule,
    RouterModule.forRoot(appRoutes,{enableTracing:true})
  ],
  providers: [ServerService,UserInputService],
  bootstrap: [AppComponent]
})
export class AppModule { }
